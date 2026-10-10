#!/usr/bin/env python3
"""Construit les données du catalogue complet à partir d'un export ColDP du Catalogue of Life.

Usage :
    python3 construire_catalogue.py <export_coldp.zip|dossier> <dossier_annotations> <dossier_sortie>

- <export_coldp.zip> : export ColDP de ChecklistBank (NameUsage.tsv, VernacularName.tsv).
- <dossier_annotations> : fichiers texte de l'atlas commenté (noms français et notes).
- <dossier_sortie> : reçoit meta.json, t/*.gz (arbre) et s/*.gz (index de recherche).

Format de sortie (texte UTF-8 compressé en gzip) :
- Les nœuds sont numérotés dans l'ordre d'un parcours en profondeur ; l'identifiant d'un
  enfant se déduit de celui du parent et de la taille des sous-arbres des aînés.
- t/N.gz : blocs « #<id du parent> » suivis d'une ligne par enfant :
    taxon :  code \t nom \t taille \t nb_espèces \t nom_fr \t note \t auteur
    espèce : ~épithète \t auteur \t nom_fr \t note   (nom complet = genre + épithète)
    espèce hors genre : !nom complet \t auteur \t nom_fr \t note
- s/<préfixe>.gz : clé normalisée \t id \t code \t nom \t nom_fr
"""
import csv, gzip, io, json, os, re, sys, unicodedata, zipfile, collections, datetime

csv.field_size_limit(1 << 30)

RANKS = {  # rang ColDP -> code de l'atlas
    "realm": "D", "domain": "D", "kingdom": "R", "subkingdom": "Rb",
    "phylum": "P", "subphylum": "Pb", "superclass": "Cs", "class": "C", "subclass": "Cb",
    "infraclass": "Ci", "superorder": "Os", "order": "O", "suborder": "Ob", "infraorder": "Oi",
    "superfamily": "Fs", "family": "F", "subfamily": "Fb", "tribe": "T", "genus": "G", "species": "S",
}
ORDER = {c: i for i, c in enumerate("D R Rb P Pb Cs C Cb Ci Os O Ob Oi Fs F Fb T G S x".split())}
ACCEPTED = {"accepted", "provisionally accepted"}

# Végétal multicellulaire exclu ; les lignées unicellulaires de Plantae sont gardées.
EXCLUDE = {
    "Tracheophyta": "Plantes vasculaires",
    "Bryophyta": "Mousses",
    "Marchantiophyta": "Hépatiques",
    "Anthocerotophyta": "Anthocérotes",
    "Florideophyceae": "Algues rouges multicellulaires (floridées)",
    "Bangiophyceae": "Algues rouges multicellulaires (nori…)",
    "Compsopogonophyceae": "Algues rouges filamenteuses",
    "Ulvophyceae": "Algues vertes multicellulaires (ulves, caulerpes…)",
    "Charophyceae": "Characées",
    "Coleochaetophyceae": "Coléochétophycées",
}
KINGDOM_FR = {
    "Animalia": "Animaux", "Fungi": "Champignons", "Plantae": "Lignée verte et rouge (unicellulaires)",
    "Chromista": "Chromistes", "Protozoa": "Protozoaires", "Bacteria": "Bactéries",
    "Archaea": "Archées", "Viruses": "Virus",
}
ALIASES = {"Metazoa": "Animalia", "Virus": "Viruses"}


def norm(s):
    s = unicodedata.normalize("NFD", s)
    s = "".join(c for c in s if unicodedata.category(c) != "Mn").lower()
    s = re.sub(r"[^a-z0-9]+", " ", s)
    return s.strip()


def open_tsv(src, wanted):
    """Ouvre NameUsage.tsv ou VernacularName.tsv dans un zip ou un dossier."""
    if os.path.isdir(src):
        for root, _, files in os.walk(src):
            for f in files:
                if f.lower() == wanted.lower():
                    return open(os.path.join(root, f), encoding="utf-8", newline="")
        return None
    z = zipfile.ZipFile(src)
    for n in z.namelist():
        if n.split("/")[-1].lower() == wanted.lower():
            return io.TextIOWrapper(z.open(n), encoding="utf-8", newline="")
    return None


def rows(fh):
    r = csv.reader(fh, delimiter="\t", quoting=csv.QUOTE_NONE)
    head = [h.split(":")[-1] for h in next(r)]
    idx = {h: i for i, h in enumerate(head)}
    for row in r:
        yield lambda k, row=row: row[idx[k]] if k in idx and idx[k] < len(row) else ""


def load_annotations(folder):
    ann = {}
    if not folder or not os.path.isdir(folder):
        return ann
    for f in sorted(os.listdir(folder)):
        if not f.endswith(".txt"):
            continue
        for line in open(os.path.join(folder, f), encoding="utf-8"):
            s = line.strip()
            if not s or s.startswith("#") or s.startswith("@@"):
                continue
            p = [x.strip() for x in s.split("|")]
            if len(p) < 2 or p[0] == "x":
                continue
            name = p[1].lstrip("†").strip()
            name = ALIASES.get(name, name)
            ann.setdefault(name, (p[2] if len(p) > 2 else "", p[3] if len(p) > 3 else ""))
    return ann


def main(src, ann_dir, out):
    ann = load_annotations(ann_dir)
    print("annotations", len(ann), file=sys.stderr)

    parent, rank, name, auth = {}, {}, {}, {}
    fh = open_tsv(src, "NameUsage.tsv")
    if fh is None:
        sys.exit("NameUsage.tsv introuvable")
    n = 0
    for g in rows(fh):
        n += 1
        if g("status").lower() not in ACCEPTED:
            continue
        i = g("ID")
        parent[i] = g("parentID")
        rank[i] = g("rank").lower()
        name[i] = g("scientificName")
        auth[i] = g("authorship")
    print("lignes", n, "acceptées", len(parent), file=sys.stderr)

    fr = {}
    vh = open_tsv(src, "VernacularName.tsv")
    if vh is not None:
        for g in rows(vh):
            if g("language").lower() in ("fra", "fre", "fr") and g("taxonID") in parent:
                fr.setdefault(g("taxonID"), g("name").strip())
    print("noms français COL", len(fr), file=sys.stderr)

    # Rattache chaque taxon gardé à son plus proche ancêtre gardé (on saute sous-genres, sections…).
    kept = {i for i, r in rank.items() if r in RANKS}
    memo = {}
    def kept_parent(i):
        chain = []
        p = parent.get(i, "")
        while p and p in parent and p not in kept:
            if p in memo:
                p = memo[p]
                break
            chain.append(p)
            p = parent.get(p, "")
        res = p if (p in kept) else ""
        for c in chain:
            memo[c] = res
        return res

    kids = collections.defaultdict(list)
    roots = []
    for i in kept:
        p = kept_parent(i)
        (kids[p].append(i) if p else roots.append(i))

    # Racine virtuelle et regroupement des eucaryotes.
    VIRT = {"__vie": ("c", "Vie", "Le vivant", "Catalogue complet des espèces décrites, d'après le Catalogue of Life. Seules les espèces acceptées sont listées (pas les synonymes ni les sous-espèces)."),
            "__euk": ("D", "Eukaryota", "Eucaryotes", "Cellules à noyau.")}
    by_name = {name[r]: r for r in roots}
    top = []
    for k in ("Bacteria", "Archaea"):
        if k in by_name: top.append(by_name[k])
    euk = [by_name[k] for k in ("Protozoa", "Chromista", "Fungi", "Animalia", "Plantae") if k in by_name]
    top.append("__euk")
    kids["__euk"] = euk
    if "Viruses" in by_name: top.append(by_name["Viruses"])
    others = [r for r in roots if r not in top and r not in euk]
    kids["__vie"] = top + others
    print("racines", [name[r] for r in roots], file=sys.stderr)

    # Exclusions végétales : sous-arbre remplacé par un nœud « x ».
    for i in list(kept):
        if name[i] in EXCLUDE and rank[i] in ("phylum", "class"):
            kids[i] = []
            rank[i] = "__x"
            fr[i] = EXCLUDE[name[i]]

    def code(i):
        if i in VIRT: return VIRT[i][0]
        if rank[i] == "__x": return "x"
        return RANKS[rank[i]]

    def sort_key(i):
        return (ORDER[code(i)], (name[i] if i not in VIRT else "").lower())

    # Parcours en profondeur itératif : identifiants, tailles, nombre d'espèces.
    order, ident = [], {}
    stack = ["__vie"]
    while stack:
        i = stack.pop()
        ident[i] = len(order)
        order.append(i)
        ch = kids.get(i, []) if i in VIRT else sorted(kids.get(i, []), key=sort_key)
        kids[i] = ch
        stack.extend(reversed(ch))
    size, nsp = {}, {}
    for i in reversed(order):
        ch = kids.get(i, [])
        size[i] = sum(1 + size[c] for c in ch)
        nsp[i] = (1 if code(i) == "S" else 0) + sum(nsp[c] for c in ch)
    print("nœuds", len(order), "espèces", nsp["__vie"], file=sys.stderr)

    def info(i):
        if i in VIRT:
            c, nm, f, note = VIRT[i]
            return c, nm, f, note, ""
        nm = name[i]
        a = ann.get(nm, ("", ""))
        f = a[0] or fr.get(i, "") or (KINGDOM_FR.get(nm, "") if rank[i] == "kingdom" else "")
        note = a[1]
        if code(i) == "x":
            note = "Végétal multicellulaire : hors du périmètre de l'atlas."
        return code(i), nm, f, note, auth.get(i, "")

    clean = lambda s: s.replace("\t", " ").replace("\n", " ")
    os.makedirs(os.path.join(out, "t"), exist_ok=True)
    os.makedirs(os.path.join(out, "s"), exist_ok=True)

    # Morceaux de l'arbre.
    chunks, buf, start = [], [], None
    def flush():
        if buf:
            k = len(chunks)
            with gzip.open(os.path.join(out, "t", f"{k}.gz"), "wt", encoding="utf-8", compresslevel=9) as w:
                w.write("".join(buf))
            chunks.append(start)
            buf.clear()
    LIMIT = 1_200_000
    cur = 0
    for i in order:
        ch = kids.get(i, [])
        if not ch:
            continue
        if start is None: start = ident[i]
        lines = [f"#{ident[i]}\n"]
        genus = name[i] + " " if code(i) == "G" and i not in VIRT else None
        for c in ch:
            cc, nm, f, note, au = info(c)
            if cc == "S":
                if genus and nm.startswith(genus):
                    lines.append(f"~{clean(nm[len(genus):])}\t{clean(au)}\t{clean(f)}\t{clean(note)}\n")
                else:
                    lines.append(f"!{clean(nm)}\t{clean(au)}\t{clean(f)}\t{clean(note)}\n")
            else:
                lines.append(f"{cc}\t{clean(nm)}\t{size[c]}\t{nsp[c]}\t{clean(f)}\t{clean(note)}\t{clean(au)}\n")
        block = "".join(lines)
        if cur + len(block) > LIMIT and buf:
            flush(); start = ident[i]; cur = 0
        buf.append(block); cur += len(block)
    flush()

    # Index de recherche : noms scientifiques et mots des noms français.
    entries = []
    for i in order:
        if i in VIRT: continue
        cc, nm, f, note, au = info(i)
        if cc == "x": continue
        k = norm(nm)
        if k: entries.append((k, ident[i], cc, nm, f))
        if f:
            words = norm(f).split()
            for w in range(len(words)):
                kk = " ".join(words[w:])
                if len(kk) >= 2: entries.append((kk, ident[i], cc, nm, f))
    entries.sort(key=lambda e: (e[0], ORDER[e[2]], e[1]))
    print("entrées de recherche", len(entries), file=sys.stderr)

    SHARD = 450_000
    def line(e): return f"{e[0]}\t{e[1]}\t{e[2]}\t{clean(e[3])}\t{clean(e[4])}\n"
    shards = {}
    def split(items, depth):
        groups = collections.OrderedDict()
        for e in items:
            p = e[0][:depth].replace(" ", "_")
            groups.setdefault(p, []).append(e)
        for p, g in groups.items():
            sz = sum(len(line(e)) for e in g)
            if sz > SHARD and depth < 6 and any(len(e[0]) > depth for e in g):
                exact = [e for e in g if len(e[0]) <= depth]
                rest = [e for e in g if len(e[0]) > depth]
                if exact: shards.setdefault(p + "=", exact)
                split(rest, depth + 1)
            else:
                shards[p] = g
    split(entries, 2)
    names = sorted(shards)
    for idx, p in enumerate(names):
        with gzip.open(os.path.join(out, "s", f"{idx}.gz"), "wt", encoding="utf-8", compresslevel=9) as w:
            w.write("".join(line(e) for e in shards[p]))

    counts = {}
    for c in kids["__vie"] + kids["__euk"]:
        if c != "__euk": counts[name[c]] = nsp[c]
    meta = {
        "built": datetime.date.today().isoformat(),
        "source": os.path.basename(src.rstrip("/")),
        "nodes": len(order), "species": nsp["__vie"], "counts": counts,
        "chunks": chunks, "shards": names,
    }
    with open(os.path.join(out, "meta.json"), "w", encoding="utf-8") as w:
        json.dump(meta, w, ensure_ascii=False, separators=(",", ":"))
    print(json.dumps({k: v for k, v in meta.items() if k not in ("chunks", "shards")}, ensure_ascii=False), file=sys.stderr)
    print("morceaux", len(chunks), "index", len(names), file=sys.stderr)


if __name__ == "__main__":
    main(*sys.argv[1:4])
