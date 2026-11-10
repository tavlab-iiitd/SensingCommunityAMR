"""2024 urinary resistance comparisons and surveillance footprints."""

from pathlib import Path
import hashlib
import json
import re
import unicodedata

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/surveillance"
DATA = OUT / "source_data"
DATA.mkdir(exist_ok=True, parents=True)
from name_cleaning import clean_organism, clean_antibiotic

YEAR = 2024
NCDC_URL = "https://ncdc.mohfw.gov.in/uploads/pdf/amr30.pdf"
ICMR_URL = "https://www.icmr.gov.in/icmrobject/uploads/Report/1763981012_icmramrsnannualreport2024.pdf"


NCDC = {
 "Escherichia coli": {
  "Ampicillin": (30771,87), "Amoxicillin/Clavulanic acid": (31408,59),
  "Piperacillin/Tazobactam": (34001,40), "Ceftriaxone": (22108,74),
  "Cefotaxime": (34490,76), "Cefepime": (28581,54), "Ertapenem": (17191,25),
  "Imipenem": (32038,23), "Meropenem": (26633,20), "Amikacin": (34617,25),
  "Gentamicin": (30890,32), "Ciprofloxacin": (36167,73),
  "Trimethoprim/Sulfamethoxazole": (37565,55), "Fosfomycin": (21602,4),
  "Nitrofurantoin": (41460,19)},
 "Klebsiella pneumoniae": {
  "Amoxicillin/Clavulanic acid": (15375,64), "Piperacillin/Tazobactam": (16141,48),
  "Ceftriaxone": (10801,68), "Cefotaxime": (15466,70), "Cefepime": (13777,56),
  "Ertapenem": (8119,41), "Imipenem": (15212,32), "Meropenem": (13196,31),
  "Amikacin": (16377,40), "Gentamicin": (15092,40), "Ciprofloxacin": (17245,61),
  "Trimethoprim/Sulfamethoxazole": (17110,52), "Nitrofurantoin": (18814,49)},
 "Pseudomonas aeruginosa": {
  "Amikacin": (5184,41), "Aztreonam": (3947,38), "Ceftazidime": (5304,55),
  "Ciprofloxacin": (4945,58), "Gentamicin": (3411,41), "Imipenem": (4949,41),
  "Meropenem": (3939,38), "Netilmicin": (3057,45), "Piperacillin/Tazobactam": (5068,32)}
}

NCDC_OPD = {
 "Ampicillin":84, "Amoxicillin/Clavulanic acid":52, "Piperacillin/Tazobactam":32,
 "Cefuroxime":68, "Ceftriaxone":66, "Cefotaxime":69, "Cefepime":43,
 "Ertapenem":17, "Imipenem":18, "Meropenem":15, "Amikacin":20,
 "Gentamicin":26, "Ciprofloxacin":69, "Trimethoprim/Sulfamethoxazole":50,
 "Fosfomycin":4, "Nitrofurantoin":17
}

ICMR = {
 "Escherichia coli": {
  "Piperacillin/Tazobactam":(5182,9269), "Cefazolin":(453,3727),
  "Cefotaxime":(1913,7422), "Ceftazidime":(1600,4722), "Ertapenem":(6975,9066),
  "Imipenem":(5989,8708), "Meropenem":(7476,10241), "Amikacin":(10221,13476),
  "Ciprofloxacin":(3020,13476), "Levofloxacin":(951,4330),
  "Trimethoprim/Sulfamethoxazole":(4486,9094), "Fosfomycin":(6288,6821),
  "Nitrofurantoin":(8981,11119)},
 "Klebsiella pneumoniae": {
  "Piperacillin/Tazobactam":(1182,3499), "Cefazolin":(169,1283),
  "Cefotaxime":(799,2764), "Ceftazidime":(468,1858), "Ertapenem":(1637,3305),
  "Imipenem":(1388,3244), "Meropenem":(1981,4078), "Amikacin":(2435,4871),
  "Ciprofloxacin":(1363,4850), "Levofloxacin":(415,1656),
  "Trimethoprim/Sulfamethoxazole":(1537,3417), "Fosfomycin":(1499,2530),
  "Nitrofurantoin":(736,3479)},
 "Pseudomonas aeruginosa": {
  "Piperacillin/Tazobactam":(667,1175), "Cefepime":(592,1140),
  "Ceftazidime":(707,1486), "Imipenem":(533,1175), "Meropenem":(609,1335),
  "Amikacin":(887,1505), "Levofloxacin":(343,921)}
}
NCDC_SITES = """BJ Medical College Ahmedabad|Gujarat
BJ Medical College Pune|Maharashtra
Government Medical College and Hospital Chandigarh|Chandigarh
GSVM Medical College Kanpur|Uttar Pradesh
Lady Hardinge Medical College and Associated Hospitals|Delhi
Mysore Medical College|Karnataka
SMS Medical College Jaipur|Rajasthan
Vardhman Mahavir Medical College and Safdarjung Hospital|Delhi
Government Medical College Thiruvananthapuram|Kerala
KAPV Government Medical College Tiruchirappalli|Tamil Nadu
Gauhati Medical College and Hospital|Assam
NEIGRIHMS Shillong|Meghalaya
MGM College and Hospital Indore|Madhya Pradesh
Indira Gandhi Medical College Shimla|Himachal Pradesh
Government Medical College and Hospital Aurangabad|Maharashtra
Osmania Medical College Hyderabad|Telangana
Guntur Medical College|Andhra Pradesh
Agartala Government Medical College|Tripura
SCB Medical College and Hospital Cuttack|Odisha
Government Medical College and Hospital Jammu|Jammu and Kashmir
PGIMS Rohtak|Haryana
Rajendra Institute of Medical Sciences Ranchi|Jharkhand
Indira Gandhi Institute of Medical Sciences Patna|Bihar
Government Medical College Haldwani|Uttarakhand
Pt Jawaharlal Nehru Memorial Medical College Raipur|Chhattisgarh
Gandhi Medical College Bhopal|Madhya Pradesh
Calcutta School of Tropical Medicine|West Bengal
GMERS Medical College and Civil Hospital Valsad|Gujarat
LLRM Medical College Meerut|Uttar Pradesh
Coimbatore Medical College and Hospital|Tamil Nadu
Maulana Azad Medical College and Associated Hospitals|Delhi
Sardar Patel Medical College and Hospital Bikaner|Rajasthan
Karnataka Institute of Medical Sciences Hubli|Karnataka
Indira Gandhi Medical College and Research Institute Puducherry|Puducherry
NAMO Medical Education and Research Institute Silvassa|Dadra and Nagar Haveli and Daman and Diu
Goa Medical College and Hospital Bambolim|Goa
STNM Medical College and Hospital Gangtok|Sikkim
Government Medical College Patiala|Punjab
Zoram Medical College Falkawn|Mizoram
Andaman and Nicobar Islands Institute of Medical Sciences|Andaman and Nicobar Islands
Jawahar Lal Nehru Institute of Medical Sciences|Manipur
Government Medical College Srinagar|Jammu and Kashmir
Rabindranath Tagore Medical College Udaipur|Rajasthan
Andhra Medical College Visakhapatnam|Andhra Pradesh
Vijayanagar Institute of Medical Sciences Ballari|Karnataka
Burdwan Medical College and Hospital|West Bengal
Grant Government Medical College and Sir JJ Group of Hospitals Mumbai|Maharashtra
Pt Raghunath Murmu Medical College and Hospital Baripada|Odisha
Government Medical College Thrissur|Kerala
SV Medical College Tirupati|Andhra Pradesh
Jorhat Medical College and Hospital|Assam
University College of Medical Sciences and GTB Hospital|Delhi
Pandit Dindayal Upadhyay Medical College Rajkot|Gujarat
Netaji Subhash Chandra Bose Medical College Jabalpur|Madhya Pradesh"""


ICMR_SITES = """CMC Vellore|Tamil Nadu
PGIMER Chandigarh|Chandigarh
JIPMER Puducherry|Puducherry
AIIMS New Delhi|Delhi
AFMC Pune|Maharashtra
AIIMS Jodhpur|Rajasthan
AIIMS Bhopal|Madhya Pradesh
Apollo Hospital Chennai|Tamil Nadu
Assam Medical College and Hospital Dibrugarh|Assam
PD Hinduja Hospital Mumbai|Maharashtra
IPGMER Kolkata|West Bengal
JPN Apex Trauma Center AIIMS Delhi|Delhi
KGMU Lucknow|Uttar Pradesh
KMC Manipal|Karnataka
LTMGH Sion Hospital Mumbai|Maharashtra
MGIMS Wardha|Maharashtra
NIMS Hyderabad|Telangana
RIMS Imphal|Manipur
SKIMS Srinagar|Jammu and Kashmir
Sir Ganga Ram Hospital New Delhi|Delhi
Tata Medical Center Kolkata|West Bengal
Artemis Hospital Gurugram|Haryana
Fortis Hospital Gurugram|Haryana"""

ALIASES = {
 "andaman nicobar islands":"andaman and nicobar islands", "chattisgarh":"chhattisgarh",
 "dadra nagar haveli":"dadra and nagar haveli and daman and diu",
 "dadra and nagar haveli":"dadra and nagar haveli and daman and diu",
 "daman diu":"dadra and nagar haveli and daman and diu",
 "daman and diu":"dadra and nagar haveli and daman and diu",
 "jammu kashmir":"jammu and kashmir", "nct delhi":"delhi", "new delhi":"delhi",
 "orissa":"odisha", "pondicherry":"puducherry", "telengana":"telangana",
 "uttaranchal":"uttarakhand"
}
def state_key(value):
    if pd.isna(value): return None
    s = unicodedata.normalize("NFKD",str(value)).encode("ascii","ignore").decode()
    s = re.sub(r"[^a-z0-9]+"," ",s.lower()).strip()
    return ALIASES.get(s,s) if s and s not in {"none","nan","unknown","not available"} else None

def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()

def prepare_data():
    roster=[]
    for network,txt,url,pages in [
        ("NCDC NARS-Net",NCDC_SITES,NCDC_URL,"57-58 (PDF 64-65)"),
        ("ICMR AMRSN",ICMR_SITES,ICMR_URL,"204 (PDF 206)")]:
        for i,line in enumerate(txt.splitlines(),1):
            institution,state=line.split("|")
            roster.append(dict(network=network,institution=institution,state=state,
                state_key=state_key(state),source_url=url,source_pages=pages,
                basis="2024 reporting laboratory" if network.startswith("NCDC") else "institution listed in 2024 report; not verified per-analysis contributor"))
    sites=pd.DataFrame(roster)
    assert sites.groupby("network").size().to_dict()=={"ICMR AMRSN":23,"NCDC NARS-Net":54}
    sites.to_csv(DATA/"network_institution_rosters_2024.csv",index=False)


    rx_path=ROOT/"data/prescriptions_clean.csv"
    rx=pd.read_csv(rx_path,usecols=["month_year","patient_rx_state","patient_rx_pincode","Rx_Count"])
    rx=rx[rx.month_year.between("2024-01","2024-12") & rx.Rx_Count.gt(0)].copy()
    rx["state_key"]=rx.patient_rx_state.map(state_key)
    dg_path=ROOT/"data/AMR_Pan_India_Aggregated_Diagnostic.csv"
    dg=pd.read_csv(dg_path)
    dg=dg[dg.end_state_date.between("2024-01-01","2024-12-31") & dg.Total_Test.gt(0)].copy()
    dg["state_key"]=dg.order_state.map(state_key)
    geo=pd.DataFrame({"prescription_records":rx.groupby("state_key").Rx_Count.sum(),
                      "diagnostic_samples":dg.groupby("state_key").Total_Test.sum()}).fillna(0)
    for network,col in [("NCDC NARS-Net","ncdc_reporting_labs"),("ICMR AMRSN","icmr_listed_institutions")]:
        geo=geo.join(sites[sites.network.eq(network)].groupby("state_key").size().rename(col),how="outer")
    geo=geo.fillna(0).astype(int).sort_index()
    geo.to_csv(DATA/"geographic_footprints_2024.csv")
    urine=dg.test_name.eq("URINE CULTURE AND SENSITIVITY")
    meta=dict(year=YEAR,prescription_records=int(rx.Rx_Count.sum()),
        prescription_states=int(rx.state_key.nunique()),prescription_raw_state_labels=int(rx.loc[rx.state_key.notna(),"patient_rx_state"].nunique()),
        prescriptions_missing_state=int(rx.loc[rx.state_key.isna(),"Rx_Count"].sum()),
        prescription_pincodes=int(rx.patient_rx_pincode.nunique()),diagnostic_samples=int(dg.Total_Test.sum()),
        diagnostic_states=int(dg.state_key.nunique()),diagnostic_pincodes=int(dg.Order_Pincode.nunique()),
        urine_samples=int(dg.loc[urine,"Total_Test"].sum()),
        diagnostic_top_five_state_percent=float(100*geo.diagnostic_samples.nlargest(5).sum()/dg.Total_Test.sum()),
        ncdc_labs=54,ncdc_states=int((geo.ncdc_reporting_labs>0).sum()),
        icmr_listed_institutions=23,icmr_roster_states=int((geo.icmr_listed_institutions>0).sum()),
        icmr_reported_analytical_centres=21,
        diagnostic_states_outside_ncdc=geo.index[(geo.diagnostic_samples>0)&(geo.ncdc_reporting_labs==0)].tolist(),
        prescription_states_outside_ncdc=geo.index[(geo.prescription_records>0)&(geo.ncdc_reporting_labs==0)].tolist(),
        prescriptions_outside_ncdc_states=int(geo.loc[geo.ncdc_reporting_labs==0,"prescription_records"].sum()))
    assert meta["prescription_records"]==998438 and meta["diagnostic_samples"]==46731
    assert meta["prescription_states"]==34 and meta["prescriptions_missing_state"]==1173
    assert meta["diagnostic_states"]==22 and meta["ncdc_states"]==32 and meta["icmr_roster_states"]==15
    (DATA/"coverage_summary_2024.json").write_text(json.dumps(meta,indent=2))

    raw_path=ROOT/"data/AMR_Pan_India_Diagnostic_Antibiotic.csv"
    raw=pd.read_csv(raw_path,usecols=["end_state_date","test_name","organism_isolated","antibiotic","RESISTANT","INTERMEDIATE","SENSITIVE","Total_Test"])
    raw=raw[raw.end_state_date.between("2024-01-01","2024-12-31") & raw.test_name.eq("URINE CULTURE AND SENSITIVITY")].copy()
    for source,dest,func in [("organism_isolated","organism_clean",clean_organism),("antibiotic","antibiotic_clean",clean_antibiotic)]:
        raw[dest]=raw[source].map({v:func(v) for v in raw[source].drop_duplicates()})
    raw=raw[raw.organism_clean.isin(ICMR)].copy()
    cnt=["RESISTANT","INTERMEDIATE","SENSITIVE"]
    assert raw[cnt].notna().all().all() and (raw[cnt]>=0).all().all()
    assert (raw[cnt].sum(axis=1)<=raw.Total_Test).all()
    ast=raw.groupby(["organism_clean","antibiotic_clean"])[cnt].sum().astype(int)
    ast["N"]=ast[cnt].sum(axis=1)
    ast["R_percent"]=100*ast.RESISTANT/ast.N
    ast["NS_percent"]=100*(ast.RESISTANT+ast.INTERMEDIATE)/ast.N


    processed_path=ROOT/"data/ast_clean.csv"
    p=pd.read_csv(processed_path,usecols=["month_year","test_name","organism_clean","antibiotic_clean",*cnt])
    p=p[p.test_name.eq("URINE") & p.organism_clean.isin(ICMR)].copy()
    pooled=p.groupby(["organism_clean","antibiotic_clean"])[cnt].sum()
    assert int(pooled.loc[("Escherichia coli","Ciprofloxacin"),"RESISTANT"])==12275
    assert int(pooled.loc[("Escherichia coli","Ciprofloxacin")].sum())==21542
    p=p[p.month_year.between("2024-01","2024-12")].groupby(["organism_clean","antibiotic_clean"])[cnt].sum().astype(int)
    pd.testing.assert_frame_equal(ast[cnt].sort_index(),p.sort_index())
    ast.to_csv(DATA/"tata_urine_ast_2024.csv")

    records=[]
    long=[]
    for org in ICMR:


        drugs=set(NCDC[org])|set(ICMR[org])
        if org=="Escherichia coli": drugs|=set(NCDC_OPD)
        for ab in sorted(drugs):
            rec=dict(year=YEAR,organism=org,antibiotic=ab)
            if (org,ab) in ast.index:
                a=ast.loc[(org,ab)]
                rec.update(tata_R=int(a.RESISTANT),tata_I=int(a.INTERMEDIATE),tata_S=int(a.SENSITIVE),
                    tata_N=int(a.N),tata_R_percent=float(a.R_percent),tata_NS_percent=float(a.NS_percent))
                for metric,numer in [("R",int(a.RESISTANT)),("NS",int(a.RESISTANT+a.INTERMEDIATE))]:
                    long.append(dict(year=YEAR,network="Tata 1mg",organism=org,antibiotic=ab,
                        specimen="urine",setting="Home-collected specimens; Tata 1mg in-house laboratories",
                        metric=metric,numerator=numer,tested_n=int(a.N),percentage=100*numer/a.N,
                        source=str(raw_path.relative_to(ROOT)),locator="2024 urinary records"))
            if ab in NCDC[org]:
                n,r=NCDC[org][ab]
                taxon={"Klebsiella pneumoniae":"Klebsiella spp.","Pseudomonas aeruginosa":"Pseudomonas spp."}.get(org,org)
                tab={"Escherichia coli":"Table 8; p.27; PDF p.34","Klebsiella pneumoniae":"Table 9; p.31; PDF p.38","Pseudomonas aeruginosa":"Table 12; p.36; PDF p.43"}[org]
                rec.update(ncdc_R_percent=r,ncdc_N=n,ncdc_taxon=taxon)
                long.append(dict(year=YEAR,network="NCDC NARS-Net",organism=taxon,antibiotic=ab,
                    specimen="urine",setting="all reported hospital locations",metric="R",numerator=None,
                    tested_n=n,percentage=r,source=NCDC_URL,locator=tab))
            if org=="Escherichia coli" and ab in NCDC_OPD:
                rec["ncdc_OPD_R_percent"]=NCDC_OPD[ab]
                long.append(dict(year=YEAR,network="NCDC NARS-Net",organism=org,antibiotic=ab,
                    specimen="urine",setting="hospital outpatient",metric="R",numerator=None,tested_n=None,
                    percentage=NCDC_OPD[ab],source=NCDC_URL,locator="Figure 12; p.29; PDF p.36; total OPD isolates 22504, not drug denominator"))
            if ab in ICMR[org]:
                s,n=ICMR[org][ab]
                assert 0<=s<=n
                rec.update(icmr_S=s,icmr_N=n,icmr_NS_percent=100*(n-s)/n)
                locator=("Table 3.2; p.100; PDF p.102" if org=="Pseudomonas aeruginosa" else "Table 2.14; p.77; PDF p.79")
                long.append(dict(year=YEAR,network="ICMR AMRSN",organism=org,antibiotic=ab,specimen="urine",
                    setting="all reported hospital locations",
                    metric="NS derived from S/N",susceptible_n=s,numerator=n-s,tested_n=n,
                    percentage=100*(n-s)/n,source=ICMR_URL,locator=locator))
            records.append(rec)
    comparison=pd.DataFrame(records)
    comparison.to_csv(DATA/"resistance_comparison_2024.csv",index=False,float_format="%.8f")
    pd.DataFrame(long).to_csv(DATA/"reported_values_and_provenance_2024.csv",index=False,float_format="%.8f")
    hashes={str(path.relative_to(ROOT)):sha256(path) for path in [raw_path,processed_path,rx_path,dg_path]}
    audit=dict(year=YEAR,source_hashes=hashes,raw_and_processed_AST_counts_match=True,
        manuscript_pooled_E_coli_ciprofloxacin_reproduced=True,
        reference_access_date="2026-10-24",n_comparison_rows=len(comparison),
        notes=["No inferential cross-network tests or pooled estimates.",
               "ICMR NS calculated from integer S/N, not from rounded percentages.",
               "NCDC percentages retained at published precision; resistant counts not reconstructed.",
               "Tata 1mg is aggregated and cannot be deduplicated at patient level."])
    (DATA/"validation_and_sources.json").write_text(json.dumps(audit,indent=2))
    return comparison,meta


if __name__ == "__main__":
    comparison, metadata = prepare_data()
    print(f"{len(comparison)} organism–antibiotic comparisons")
