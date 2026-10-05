"""Name standardisation shared by the preprocessing scripts.

1. WHO AWaRe 2025 matching of prescription products and laboratory agents.
2. Cleaning of free-text organism and antibiotic names in the laboratory data.
"""

import difflib
import re

import pandas as pd

WHO_FILE = "data/aware-2025.xlsx"

# ---------------------------------------------------------------------------
# WHO AWaRe 2025 matching
# ---------------------------------------------------------------------------

COMPONENT_SYNONYMS = {
    "amoxycillin": "amoxicillin",
    "tazobactum": "tazobactam",
    "cephalexin": "cefalexin",
    "potassium clavulanate": "clavulanic acid",
    "clavulanate": "clavulanic acid",
    "sultamicillin tosilate": "sultamicillin",
    "flucloxacillin sodium": "flucloxacillin",
    "cefditoren": "cefditoren pivoxil",
    "cefetamet": "cefetamet pivoxil",
    "benzathine penicillin g": "benzathine benzylpenicillin",
    "penicillin v": "phenoxymethylpenicillin",
    "penicillin g": "benzylpenicillin",
    "co trimoxazole": "sulfamethoxazole/trimethoprim",
    "cotrimoxazole": "sulfamethoxazole/trimethoprim",
    "serratiopeptidase": "serrapeptase",
    "bromelain": "bromelains",
    "lactobacillus sporogenes": "bacillus coagulans",
    "tobramicin": "tobramycin",          # misspelt on the WHO Not recommended sheet
    "benzyl penicillin": "benzylpenicillin",
}

WHO_LABEL_FIXES = {"tobramicin": "tobramycin", "benzyl penicillin": "benzylpenicillin"}

# A laboratory panel tests an agent, so these panel names link to the WHO entry of that agent
AST_AGENT_LINKS = {
    "imipenem": "imipenem/cilastatin",
    "high-level gentamicin": "gentamicin",
    "high level gentamicin": "gentamicin",
}
AST_NOT_CLASSIFIED = {"penicillin", "imipenem/relebactam", "imipenem/cilastatin/relebactam"}


def components(name):
    text = re.sub(r"\s*\(.*?\)\s*", " ", str(name).strip().lower()).replace("-", " ")
    parts = []
    for part in re.split(r"[+/]", text):
        part = re.sub(r"\s+", " ", part.strip())
        if part:
            part = COMPONENT_SYNONYMS.get(part, part)
            parts += [p.strip() for p in part.split("/") if p.strip()]
    return parts


def match_key(name):
    """Order-insensitive key, so A+B and B+A match the same WHO entry."""
    return "/".join(sorted(components(name)))


def load_who(path=WHO_FILE):
    awr = pd.read_excel(path, sheet_name="AWaRe classification 2025", header=3)
    awr = awr.dropna(subset=["Antibiotic"])
    sheet = pd.read_excel(path, sheet_name="Not recommended", header=None)
    not_recommended = [str(v).strip() for v in sheet[0]]
    not_recommended = [v for v in not_recommended
                       if v not in ("nan", "", "Antibiotic") and len(v) < 90
                       and not v.lower().startswith(("the use of", "antibiotics not recommended"))]

    category, who_class = {}, {}
    for _, row in awr.iterrows():
        name = str(row.Antibiotic).strip()
        key = match_key(re.sub(r"_(IV|oral)$", "", name, flags=re.I))
        if name.lower().endswith("_oral") or key not in category:   # oral route wins
            category[key] = str(row.Category).strip()
            who_class[key] = str(row.Class).strip().replace("_anti-pseudomonal", " (anti-pseudomonal)")

    label = {}
    for name in not_recommended + [re.sub(r"_(IV|oral)$", "", n.strip(), flags=re.I)
                                   for n in awr["Antibiotic"].astype(str)]:
        label.setdefault(match_key(name), name)

    return {"category": category, "class": who_class, "label": label,
            "not_recommended": {match_key(n) for n in not_recommended}}


def aware_group(name, who):
    """(AWaRe group, WHO class) for a whole product: Not recommended first, then A/W/R."""
    key = match_key(name)
    if key in who["not_recommended"]:
        return "Not recommended", "Not recommended FDC"
    if key in who["category"]:
        return who["category"][key], who["class"][key]
    return "Unclassified", None


def lab_aware_group(agent, who):
    key = str(agent).strip().lower()
    if key in AST_NOT_CLASSIFIED:
        return "Unclassified", None
    return aware_group(AST_AGENT_LINKS.get(key, agent), who)


def product_label(product, who):
    """WHO spelling for matched products, INN spelling for unmatched ones."""
    name = who["label"].get(match_key(product))
    if name is not None:
        parts = [re.sub(r"\s+", " ", p.strip()) for p in name.split("/")]
        name = "/".join(WHO_LABEL_FIXES.get(p.lower(), p) for p in parts)
    else:
        parts = []
        for part in re.split(r"[+/]", str(product)):
            text = re.sub(r"\s+", " ", part.strip())
            if text:
                key = re.sub(r"\s*\(.*?\)\s*", " ", text.lower()).replace("-", " ")
                parts.append(COMPONENT_SYNONYMS.get(re.sub(r"\s+", " ", key).strip(), text.lower()))
        name = "/".join(parts)
    return name[0].upper() + name[1:] if name else name


# ---------------------------------------------------------------------------
# Organism and antibiotic names in the laboratory data
# ---------------------------------------------------------------------------

CANONICAL_ORGANISMS = [
    "Escherichia coli",
    "Klebsiella pneumoniae", "Klebsiella oxytoca", "Klebsiella species",
    "Acinetobacter baumannii", "Acinetobacter baumannii complex", "Acinetobacter species",
    "Pseudomonas aeruginosa", "Pseudomonas species",
    "Enterobacter cloacae", "Enterobacter cloacae complex", "Enterobacter species",
    "Citrobacter freundii", "Citrobacter koseri", "Citrobacter species",
    "Proteus mirabilis", "Proteus vulgaris", "Proteus species",
    "Morganella morganii",
    "Serratia marcescens", "Serratia species",
    "Stenotrophomonas maltophilia",
    "Haemophilus influenzae",
    "Salmonella Typhi", "Salmonella Paratyphi", "Salmonella species",
    "Shigella sonnei", "Shigella species",
    "Providencia rettgeri", "Providencia stuartii", "Providencia species",
    "Aeromonas hydrophila", "Aeromonas species",
    "Burkholderia cepacia",
    "Staphylococcus aureus",
    "Coagulase negative Staphylococcus",
    "Staphylococcus species",
    "Streptococcus pneumoniae", "Streptococcus pyogenes",
    "Streptococcus agalactiae", "Streptococcus species",
    "Enterococcus faecalis", "Enterococcus faecium", "Enterococcus species",
    "Candida albicans", "Candida parapsilosis", "Candida tropicalis", "Candida species",
    "Vibrio alginolyticus",
    "Pluralibacter gergoviae",
    "Achromobacter xylosoxidans",
    "Sphingomonas paucimobilis",
    "Lelliottia amnigena",
    "Corynebacterium diphtheriae",
]

CANONICAL_ANTIBIOTICS = [
    "Amikacin", "Amoxicillin", "Amoxicillin/Clavulanic acid",
    "Amphotericin B", "Ampicillin", "Ampicillin/Sulbactam",
    "Azithromycin", "Aztreonam",
    "Benzylpenicillin",
    "Caspofungin",
    "Cefazolin", "Cefepime", "Cefixime", "Cefoperazone", "Cefoperazone/Sulbactam",
    "Cefotaxime", "Cefoxitin", "Ceftazidime", "Ceftriaxone", "Cefuroxime",
    "Chloramphenicol", "Ciprofloxacin", "Clarithromycin", "Clindamycin",
    "Colistin", "Daptomycin", "Doripenem", "Doxycycline",
    "Erythromycin", "Ertapenem",
    "Fluconazole", "Flucytosine", "Fosfomycin",
    "Gentamicin", "High-level gentamicin", "Imipenem",
    "Levofloxacin", "Linezolid",
    "Meropenem", "Micafungin", "Minocycline", "Moxifloxacin",
    "Nalidixic acid", "Netilmicin", "Nitrofurantoin", "Norfloxacin",
    "Ofloxacin", "Oxacillin", "Penicillin", "Piperacillin",
    "Piperacillin/Tazobactam", "Polymyxin B",
    "Teicoplanin", "Tetracycline",
    "Ticarcillin", "Ticarcillin/Clavulanic acid", "Tigecycline", "Tobramycin",
    "Trimethoprim/Sulfamethoxazole", "Vancomycin", "Voriconazole",
]

# Prefixes of the letters-only, upper-case name; the longest matching prefix wins
ORGANISM_PREFIXES = {
    "ESCHERICHIACOLI": "Escherichia coli",
    "KLEBSIELLAPNEUMON": "Klebsiella pneumoniae",
    "KLEBSIELLAOXYTOC": "Klebsiella oxytoca",
    "KLEBSIELLA": "Klebsiella species",
    "ACINETOBACTERBAUMANNIIC": "Acinetobacter baumannii complex",
    "ACINETOBACTERBAUMANN": "Acinetobacter baumannii",
    "ACINETOBACTER": "Acinetobacter species",
    "PSEUDOMONASAERUG": "Pseudomonas aeruginosa",
    "PSEUDOMONAS": "Pseudomonas species",
    "ENTEROBACTERCLOACAECOM": "Enterobacter cloacae complex",
    "ENTEROBACTERCLOACAE": "Enterobacter cloacae",
    "ENTEROBACTER": "Enterobacter species",
    "CITROBACTERFREUND": "Citrobacter freundii",
    "CITROBACTERKOSERI": "Citrobacter koseri",
    "CITROBACTER": "Citrobacter species",
    "PROTEUSMIRAB": "Proteus mirabilis",
    "PROTEUSVULG": "Proteus vulgaris",
    "PROTEUS": "Proteus species",
    "MORGANELLAMORGAN": "Morganella morganii",
    "SERRATIAMARCESC": "Serratia marcescens",
    "SERRATIA": "Serratia species",
    "STENOTROPHOMONASMALTOPH": "Stenotrophomonas maltophilia",
    "SALMONELLATYPHI": "Salmonella Typhi",
    "SALMONELLAPARATYP": "Salmonella Paratyphi",
    "SALMONELLA": "Salmonella species",
    "SHIGELLASONN": "Shigella sonnei",
    "SHIGELLA": "Shigella species",
    "PROVIDENCIARETTG": "Providencia rettgeri",
    "PROVIDENCIASTUAR": "Providencia stuartii",
    "PROVIDENCIA": "Providencia species",
    "AEROMONASHYDROPH": "Aeromonas hydrophila",
    "AEROMONAS": "Aeromonas species",
    "BURKHOLDERIACEP": "Burkholderia cepacia",
    "STAPHYLOCOCCUSAUR": "Staphylococcus aureus",
    "STAPHYLOCOCCUS": "Staphylococcus species",
    "STREPTOCOCCUSPNEUMON": "Streptococcus pneumoniae",
    "STREPTOCOCCUSPYOG": "Streptococcus pyogenes",
    "STREPTOCOCCUSAGAL": "Streptococcus agalactiae",
    "STREPTOCOCCUS": "Streptococcus species",
    "ENTEROCOCCUSFAECALIS": "Enterococcus faecalis",
    "ENTEROCOCCUSFAECIUM": "Enterococcus faecium",
    "ENTEROCOCCUS": "Enterococcus species",
    "CANDIDAALBIC": "Candida albicans",
    "CANDIDAPARAP": "Candida parapsilosis",
    "CANDIDATROPIC": "Candida tropicalis",
    "CANDIDA": "Candida species",
}

ANTIBIOTIC_PREFIXES = {
    "ALIDIXI": "Nalidixic acid",
    "AMIKACI": "Amikacin",
    "AMOXICIL": "Amoxicillin",
    "AMOXYCIL": "Amoxicillin",
    "AMPHOTERI": "Amphotericin B",
    "AMPICIL": "Ampicillin",
    "AZITHRO": "Azithromycin",
    "AZTREO": "Aztreonam",
    "ATREONA": "Aztreonam",
    "BENZYLPEN": "Benzylpenicillin",
    "CASPOFU": "Caspofungin",
    "CEFAZOL": "Cefazolin",
    "CEFEPI": "Cefepime",
    "CEFIPI": "Cefepime",
    "CEFEXIM": "Cefixime",
    "CEFIXIM": "Cefixime",
    "CEFOPER": "Cefoperazone",
    "CEFOXIT": "Cefoxitin",
    "CEFOTAX": "Cefotaxime",
    "CEFTAZI": "Ceftazidime",
    "CEFTRIA": "Ceftriaxone",
    "CEFUROX": "Cefuroxime",
    "CHLORAMP": "Chloramphenicol",
    "CIPROFL": "Ciprofloxacin",
    "CLARITH": "Clarithromycin",
    "CLINDAM": "Clindamycin",
    "COLISTI": "Colistin",
    "DAPTOMY": "Daptomycin",
    "DORIPE": "Doripenem",
    "DOXYCYC": "Doxycycline",
    "ERTAPEN": "Ertapenem",
    "ERYTHRO": "Erythromycin",
    "FLUCON": "Fluconazole",
    "FLUCYTO": "Flucytosine",
    "FOSFOMY": "Fosfomycin",
    "FOSOMY": "Fosfomycin",
    "GENTAMI": "Gentamicin",
    "IMIPENE": "Imipenem",
    "LEVOFLO": "Levofloxacin",
    "LINEZOL": "Linezolid",
    "MEROPEN": "Meropenem",
    "MICAFUN": "Micafungin",
    "MINOCY": "Minocycline",
    "MOXIFLO": "Moxifloxacin",
    "NALIDIX": "Nalidixic acid",
    "NETILMI": "Netilmicin",
    "NITROFU": "Nitrofurantoin",
    "NORFLOX": "Norfloxacin",
    "OFLOXAC": "Ofloxacin",
    "OXACILL": "Oxacillin",
    "PENICIL": "Penicillin",
    "PIPERAC": "Piperacillin",
    "POLYMYX": "Polymyxin B",
    "TEICOPL": "Teicoplanin",
    "TETRACY": "Tetracycline",
    "TICARCI": "Ticarcillin",
    "TIGECYC": "Tigecycline",
    "TOBRAMY": "Tobramycin",
    "TRIMETH": "Trimethoprim/Sulfamethoxazole",
    "SULFAME": "Trimethoprim/Sulfamethoxazole",
    "VANCOMY": "Vancomycin",
    "VORICON": "Voriconazole",
}

# Misspellings that the fuzzy match would send to cefixime instead of cefepime
ANTIBIOTIC_TYPOS = {"CEPEFIME": "Cefepime", "CEFEPIIME": "Cefepime",
                    "CEPEPIME": "Cefepime", "CEFEPIM": "Cefepime"}

NOVEL_BLI_PARTNERS = ["RELEBACTAM", "VABORBACTAM", "AVIBACTAM", "TANIBORBACTAM", "ZIDEBACTAM"]

ORGANISM_TYPOS = {
    r'acinetrobacter': 'acinetobacter',
    r'enterbacter': 'enterobacter',
    r'entrobacter': 'enterobacter',
    r'klebis?ella': 'klebsiella',
    r'citrobater|citribacter': 'citrobacter',
    r'morganella\s+morgani\b': 'morganella morganii',
    r'aeriginosa': 'aeruginosa',
    r'escherichiacoli': 'escherichia coli',
    r'\bpeumoniae?\b': 'pneumoniae',
    r'\bpnuemoniae?\b': 'pneumoniae',
    r'\bpnemoniae?\b': 'pneumoniae',
    r'\bpnemonoae?\b': 'pneumoniae',
    r'\bcolacae\b': 'cloacae',
    r'\bvulagris\b': 'vulgaris',
}


def letters(text):
    return re.sub(r'[^A-Z]', '', text.upper())


def prefix_lookup(text, prefixes):
    matches = [(k, v) for k, v in prefixes.items() if text.startswith(k)]
    return max(matches, key=lambda kv: len(kv[0]))[1] if matches else None


def closest(name, choices):
    match = difflib.get_close_matches(name, choices, n=1, cutoff=0.75)
    return match[0] if match else None


def clean_organism(raw):
    if pd.isna(raw) or str(raw).strip() == "":
        return pd.NA

    s = str(raw).lower().strip()
    s = re.split(r'antibiogram|mic interpretation|antibiotic susceptibility', s)[0]
    s = re.sub(r'note\s*:?|kindly correlate clinically|comment confirmed'
               r'|\bisolated\b|\bseen\b|\band\s*$', '', s)
    s = re.sub(r'\s+', ' ', s).strip()
    if not s or re.search(r'^\d|\bcfu\b', s):           # colony counts, not organisms
        return pd.NA

    if re.fullmatch(r'no\s*(significant\s*)?growth', s) or \
       re.search(r'\bno\s*pathogen(ic)?\s*organism\b', s):
        return "No growth"

    for pattern, correction in ORGANISM_TYPOS.items():
        s = re.sub(pattern, correction, s)

    if 'staphylococcus aureus' in s or re.search(r'\bmrsa\b|\bmssa\b', s):
        return "Staphylococcus aureus"
    if ('coagulase' in s and 'negative' in s) or re.search(r'\bmr-?cons\b', s):
        return "Coagulase negative Staphylococcus"
    if 'enterococcus' in s and 'vancomycin' in s:
        return "Enterococcus species"
    if re.fullmatch(r'methicill?in(\s+resistant)?', s):
        return pd.NA
    if 'enterobacter' in s and 'cloacae' in s and 'complex' in s:
        return "Enterobacter cloacae complex"
    if 'acinetobacter' in s and 'baumannii' in s and 'complex' in s:
        return "Acinetobacter baumannii complex"

    s = re.sub(r'\bsps\b|\bspp\.?\b|\bsp\.\b', 'species', s)
    s = re.sub(r'\s+', ' ', s).strip()

    tokens = s.split()
    base = ' '.join(tokens[:2]) if len(tokens) >= 2 else s     # genus + species only
    return (prefix_lookup(letters(base), ORGANISM_PREFIXES)
            or closest(base.title(), CANONICAL_ORGANISMS)
            or closest(s.title(), CANONICAL_ORGANISMS)
            or s[0].upper() + s[1:])


def clean_antibiotic(raw):
    if pd.isna(raw) or str(raw).strip() == "":
        return pd.NA

    raw = str(raw).strip()
    lab_code = re.match(r'^[a-z]+([A-Z].*)$', raw)              # e.g. "azOFLOXACIN"
    s = lab_code.group(1).upper() if lab_code else raw.upper()
    s = re.sub(r'(?<=[A-Z])0(?=[A-Z])', 'O', s)                 # AZTRE0NAM
    s = re.sub(r'(?<=[A-Z])6(?=[A-Z])', 'G', s)
    s = re.sub(r'[\d\\*,`]', ' ', s)
    s = re.sub(r'\b(SENSITIVE|RESISTANT|INTERMEDIATE)\b', '', s)
    s = re.sub(r'^(SENSITIVE|RESISTANT|INTERMEDIATE)\s*', '', s)
    s = re.sub(r'\s+', ' ', s).strip()
    s_letters = letters(s)

    for token in NOVEL_BLI_PARTNERS:                            # keep e.g. imipenem/relebactam separate
        if token in s_letters:
            base = prefix_lookup(s_letters, ANTIBIOTIC_PREFIXES)
            return f"{base}/{token.lower()}" if base else s.title()

    if re.search(r'HIGH\s*LEVEL', s) and 'GENT' in s:
        return "High-level gentamicin"
    if 'AMOX' in s and 'CLAV' in s:
        return "Amoxicillin/Clavulanic acid"
    if 'PIPERA' in s and 'TAZO' in s:
        return "Piperacillin/Tazobactam"
    if 'CEFOPER' in s and 'SULB' in s:
        return "Cefoperazone/Sulbactam"
    if 'TICAR' in s and 'CLAV' in s:
        return "Ticarcillin/Clavulanic acid"
    if 'AMPICILLIN' in s and 'SULBACT' in s:
        return "Ampicillin/Sulbactam"
    if re.search(r'\bCO-?TRIM', s) or ('TRIMETHOPRIM' in s and 'SULFA' in s):
        return "Trimethoprim/Sulfamethoxazole"

    return (prefix_lookup(s_letters, ANTIBIOTIC_PREFIXES)
            or ANTIBIOTIC_TYPOS.get(s_letters)
            or closest(s.title(), CANONICAL_ANTIBIOTICS)
            or s.title())
