#!/usr/bin/env python3
"""Construit les données du catalogue complet à partir d'un export ColDP du Catalogue of Life.

Usage :
    python3 construire_catalogue.py <export_coldp.zip|dossier> <dossier_annotations> <dossier_sortie> [<version_de_base.zip>]

- <export_coldp.zip> : export ColDP de ChecklistBank (NameUsage.tsv, VernacularName.tsv).
- <dossier_annotations> : fichiers texte de l'atlas commenté (noms français et notes).
- <dossier_sortie> : reçoit meta.json, t/*.gz (arbre) et s/*.gz (index de recherche).
- <version_de_base.zip> (facultatif) : export ColDP de la version de base, d'où sont repris noms
  français et statut « éteint » quand l'export principal (version étendue) ne les contient pas.

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
    "Pteridobiotina": "Plantes vasculaires (fougères, conifères, plantes à fleurs…)",
    "Bryobiotina": "Mousses, hépatiques et anthocérotes, avec leurs fossiles",
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


def usages(src):
    """Produit (id, parent, statut, rang, nom, auteur) quel que soit le format d'export."""
    fh = open_tsv(src, "NameUsage.tsv")
    if fh is not None:
        print("format : ColDP NameUsage", file=sys.stderr)
        for g in rows(fh):
            yield g("ID"), g("parentID"), g("status") or "accepted", g("rank"), g("scientificName"), g("authorship"), g("extinct")
        return
    fh = open_tsv(src, "Taxon.tsv")
    if fh is None:
        sys.exit("Ni NameUsage.tsv ni Taxon.tsv dans l'export")
    it = rows(fh)
    first = next(it, None)
    if first is None:
        return
    if first("taxonID"):
        print("format : Darwin Core", file=sys.stderr)
        for g in _chain(first, it):
            yield (g("taxonID"), g("parentNameUsageID"), g("taxonomicStatus") or "accepted", g("taxonRank"),
                   g("scientificName"), g("scientificNameAuthorship"), "")
        return
    print("format : ColDP Taxon + Name", file=sys.stderr)
    taxa = []
    for g in _chain(first, it):
        st = "provisionally accepted" if g("provisional").lower() == "true" else "accepted"
        taxa.append((g("ID"), g("parentID"), st, g("nameID"), g("extinct")))
    need = {t[3] for t in taxa}
    names = {}
    nh = open_tsv(src, "Name.tsv")
    for g in rows(nh):
        if g("ID") in need:
            names[g("ID")] = (g("rank"), g("scientificName"), g("authorship"))
    for i, p, st, nid, ext in taxa:
        rk, nm, au = names.get(nid, ("", "", ""))
        yield i, p, st, rk, nm, au, ext


def _chain(first, it):
    yield first
    yield from it


def clean_name(nm):
    nm = nm.lstrip("†").strip()
    return ALIASES.get(nm, nm)


def load_atlas(folder):
    """Noms français et notes de l'atlas, plus ses espèces avec leurs ancêtres."""
    ann, especes, ancetres = {}, [], {}
    if not folder or not os.path.isdir(folder):
        return ann, especes
    for f in sorted(os.listdir(folder)):
        if not f.endswith(".txt"):
            continue
        base, pile = [], []
        for brut in open(os.path.join(folder, f), encoding="utf-8"):
            brut = brut.rstrip("\n")
            if not brut.strip() or brut.strip().startswith("#"):
                continue
            if brut.startswith("@@"):
                a = brut[2:].strip()
                base, pile = ancetres.get(a, []) + [a], []
                continue
            prof = len(brut) - len(brut.lstrip(" "))
            p = [x.strip() for x in brut.strip().split("|")]
            if len(p) < 2:
                continue
            code, nm = p[0], p[1]
            pile = pile[:prof]
            anc = base + pile
            ancetres.setdefault(nm, anc)
            pile.append(nm)
            if code == "x":
                continue
            ann.setdefault(clean_name(nm), (p[2] if len(p) > 2 else "", p[3] if len(p) > 3 else ""))
            if code == "S":
                especes.append((clean_name(nm), nm.startswith("†"), [clean_name(a) for a in anc]))
    return ann, especes


def main(src, ann_dir, out, base=None):
    ann, atlas_especes = load_atlas(ann_dir)
    print("annotations", len(ann), "espèces de l'atlas", len(atlas_especes), file=sys.stderr)

    parent, rank, name, auth, extinct = {}, {}, {}, {}, set()
    seen_ranks = collections.Counter()
    n = 0
    for i, p, st, rk, nm, au, ext in usages(src):
        n += 1
        if st.lower() not in ACCEPTED:
            continue
        if au and nm.endswith(" " + au):
            nm = nm[: -len(au) - 1]
        parent[i] = p
        rank[i] = rk.lower()
        name[i] = nm
        auth[i] = au
        if ext.lower() in ("true", "1", "yes"):
            extinct.add(i)
        seen_ranks[rank[i]] += 1
    print("lignes", n, "acceptées", len(parent), "éteints", len(extinct), file=sys.stderr)
    print("rangs", seen_ranks.most_common(40), file=sys.stderr)

    fr = {}
    vh = open_tsv(src, "VernacularName.tsv")
    if vh is not None:
        for g in rows(vh):
            lang = g("language").lower()
            tid = g("taxonID")
            if lang in ("fra", "fre", "fr") and tid in parent:
                fr.setdefault(tid, (g("name") or g("vernacularName")).strip())
    print("noms français COL", len(fr), file=sys.stderr)

    if base:
        vern = {}
        vh = open_tsv(base, "VernacularName.tsv")
        if vh is not None:
            for g in rows(vh):
                if g("language").lower() in ("fra", "fre", "fr"):
                    vern.setdefault(g("taxonID"), (g("name") or g("vernacularName")).strip())
        fr_nom, ext_nom = {}, set()
        for i, p, st, rk, nm, au, ext in usages(base):
            if st.lower() not in ACCEPTED:
                continue
            if au and nm.endswith(" " + au):
                nm = nm[: -len(au) - 1]
            if i in vern:
                fr_nom.setdefault(nm, vern[i])
            if ext.lower() in ("true", "1", "yes"):
                ext_nom.add(nm)
        for i, nm in name.items():
            if nm in fr_nom:
                fr.setdefault(i, fr_nom[nm])
            if nm in ext_nom:
                extinct.add(i)
        print("depuis la version de base : noms français", len(fr), "éteints", len(extinct), file=sys.stderr)
    atlas_eteints = {nm for nm, ex, _ in atlas_especes if ex}
    for i, nm in name.items():
        if nm in atlas_eteints:
            extinct.add(i)

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
    kpar = {}
    for i in kept:
        p = kept_parent(i)
        kpar[i] = p
        (kids[p].append(i) if p else roots.append(i))

    # Espèces de l'atlas absentes du catalogue : rattachées à leur genre, ou au groupe le plus proche.
    par_nom = collections.defaultdict(list)
    for i in kept:
        par_nom[name[i]].append(i)
    def anc_noms(i):
        out_ = set()
        while i:
            out_.add(name[i]); i = kpar.get(i, "")
        return out_
    def meilleur(cands, anc):
        best, score = None, 0
        for c in cands:
            sc = len(anc_noms(c) & anc)
            if sc > score:
                best, score = c, sc
        return best
    ajouts = 0
    for nm, ex, anc in atlas_especes:
        if nm in par_nom or "Viruses" in anc or " " not in nm:
            continue
        ancs = set(anc)
        genre = nm.split()[0]
        g = meilleur([c for c in par_nom.get(genre, []) if rank[c] == "genus"], ancs)
        if g is None:
            attache = None
            for a in reversed(anc):
                # Rattachement à un embranchement ou plus bas, jamais directement sous un règne.
                attache = meilleur([c for c in par_nom.get(a, []) if ORDER[RANKS[rank[c]]] >= ORDER["P"]], ancs)
                if attache:
                    break
            if attache is None:
                continue
            g = f"__g:{genre}:{attache}"
            if g not in rank:
                parent[g], rank[g], name[g], auth[g], kpar[g] = attache, "genus", genre, "", attache
                kids[attache].append(g); par_nom[genre].append(g)
        s = f"__s:{nm}"
        parent[s], rank[s], name[s], auth[s], kpar[s] = g, "species", nm, "", g
        if ex:
            extinct.add(s)
        kids[g].append(s); par_nom[nm].append(s)
        ajouts += 1
    print("espèces ajoutées depuis l'atlas", ajouts, file=sys.stderr)

    # Racine virtuelle et regroupement des eucaryotes.
    VIRT = {"__vie": ("c", "Vie", "Le vivant", "Catalogue complet des espèces décrites, d'après le Catalogue of Life. Seules les espèces acceptées sont listées (pas les synonymes ni les sous-espèces)."),
            "__euk": ("D", "Eukaryota", "Eucaryotes", "Cellules à noyau."),
            "__vir": ("D", "Viruses", "Virus (hors arbre)", "Les virus ne sont pas des cellules et n'ont pas d'ancêtre commun unique : ils sont classés à part, en « domaines » (realms) définis par l'ICTV, avec les viroïdes et les satellites.")}
    print("racines", sorted(name[r] for r in roots), file=sys.stderr)
    by_name = {name[r]: r for r in roots}
    top = [by_name[k] for k in ("Bacteria", "Archaea") if k in by_name]
    if "Eukaryota" in by_name:
        euk_node = by_name["Eukaryota"]
    else:
        euk_node = "__euk"
        kids["__euk"] = [by_name[k] for k in ("Protozoa", "Chromista", "Fungi", "Animalia", "Plantae") if k in by_name]
    top.append(euk_node)
    placed = set(top) | set(kids.get("__euk", []))
    others = [r for r in roots if r not in placed]
    # Racines sans classification : on garde les virus (anciens noms de familles…), on écarte le reste
    # (souvent des variantes orthographiques intégrées automatiquement par la version étendue).
    VIRUS_RE = re.compile(r"vir|viroid|satellit", re.I)
    ecartes = [r for r in others if not (rank[r] == "realm" or VIRUS_RE.search(name[r]))]
    others = [r for r in others if r not in ecartes]
    print("racines sans classification écartées", len(ecartes), file=sys.stderr)
    if "Viruses" in by_name:
        vir_node = by_name["Viruses"]
        kids[vir_node] = kids.get(vir_node, []) + [r for r in others if r != vir_node]
    else:
        vir_node = "__vir"
        kids["__vir"] = others
    top.append(vir_node)
    kids["__vie"] = top

    # Exclusions végétales : sous-arbre remplacé par un nœud « x ».
    for i in list(kept):
        if name[i] in EXCLUDE and rank[i] in ("subkingdom", "phylum", "class"):
            kids[i] = []
            rank[i] = "__x"
            fr[i] = EXCLUDE[name[i]]

    def code(i):
        if i in VIRT: return VIRT[i][0]
        if rank[i] == "__x": return "x"
        return RANKS[rank[i]]

    def sort_key(i):
        return (ORDER[code(i)], (VIRT[i][1] if i in VIRT else name[i]).lower())

    # Parcours en profondeur itératif : identifiants, tailles, nombre d'espèces.
    order, ident = [], {}
    stack = ["__vie"]
    while stack:
        i = stack.pop()
        ident[i] = len(order)
        order.append(i)
        ch = kids.get(i, []) if i == "__vie" else sorted(kids.get(i, []), key=sort_key)
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
            return c, nm, f, note, "", ""
        nm = name[i]
        a = ann.get(nm, ("", ""))
        f = a[0] or fr.get(i, "") or (KINGDOM_FR.get(nm, "") if rank[i] == "kingdom" else "")
        note = a[1]
        if code(i) == "x":
            note = "Végétal multicellulaire : hors du périmètre de l'atlas."
        return code(i), nm, f, note, auth.get(i, ""), ("1" if i in extinct else "")

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
            cc, nm, f, note, au, ex = info(c)
            if cc == "S":
                if genus and nm.startswith(genus):
                    lines.append(f"~{clean(nm[len(genus):])}\t{clean(au)}\t{clean(f)}\t{clean(note)}\t{ex}".rstrip("\t") + "\n")
                else:
                    lines.append(f"!{clean(nm)}\t{clean(au)}\t{clean(f)}\t{clean(note)}\t{ex}".rstrip("\t") + "\n")
            else:
                lines.append(f"{cc}\t{clean(nm)}\t{size[c]}\t{nsp[c]}\t{clean(f)}\t{clean(note)}\t{clean(au)}\t{ex}".rstrip("\t") + "\n")
        block = "".join(lines)
        if cur + len(block) > LIMIT and buf:
            flush(); start = ident[i]; cur = 0
        buf.append(block); cur += len(block)
    flush()

    # Index de recherche : noms scientifiques et mots des noms français.
    entries = []
    for i in order:
        if i in VIRT: continue
        cc, nm, f, note, au, ex = info(i)
        if cc == "x": continue
        k = norm(nm)
        if k: entries.append((k, ident[i], cc, nm, f, ex))
        if f:
            words = norm(f).split()
            for w in range(len(words)):
                kk = " ".join(words[w:])
                if len(kk) >= 2: entries.append((kk, ident[i], cc, nm, f, ex))
    entries.sort(key=lambda e: (e[0], ORDER[e[2]], e[1]))
    print("entrées de recherche", len(entries), file=sys.stderr)

    # Tranches contiguës d'environ 1 Mo ; une même clé n'est jamais coupée en deux.
    SHARD = 1_000_000
    def line(e):
        nm = "" if e[3] == e[0][:1].upper() + e[0][1:] else clean(e[3])  # nom déductible de la clé
        return f"{e[0]}\t{e[1]}\t{e[2]}\t{nm}\t{clean(e[4])}\t{e[5]}".rstrip("\t") + "\n"
    names, buf, cur, last = [], [], 0, None
    def flush_shard():
        with gzip.open(os.path.join(out, "s", f"{len(names) - 1}.gz"), "wt", encoding="utf-8", compresslevel=9) as w:
            w.write("".join(buf))
    for e in entries:
        l = line(e)
        if not buf or (cur + len(l) > SHARD and e[0] != last):
            if buf: flush_shard()
            names.append(e[0]); buf, cur = [], 0
        buf.append(l); cur += len(l); last = e[0]
    if buf: flush_shard()

    counts = {}
    for c in kids["__vie"] + kids.get(euk_node, []):
        if c != euk_node: counts[c if c in VIRT else name[c]] = nsp[c]
    counts = {(VIRT[k][1] if k in VIRT else k): v for k, v in counts.items()}
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
    try:
        import resource
        print("mémoire max (Mo)", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024, file=sys.stderr)
    except ImportError:
        pass


if __name__ == "__main__":
    main(*sys.argv[1:5])
