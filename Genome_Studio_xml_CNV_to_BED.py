# Tool to combine all xml files from directory to single BED file.
# Input: Directory with xml files, patient data from excel file, text files to keep track of sample numbers
# Output: xlsx-file and BED-fail of all automatically detected CNV-s.
# 16.01.2024
# v4 (xml based)
# Kadi Jairus


import os
import pandas as pd
from datetime import datetime

print("Programm alustas tööd, palun oota...")

# Serveris asuvate asjade asukohad
# koondtabeli_asukoht = r'\\srvlaste\Yhendlabor\GE_Illumina kiip\GDA_38_koond.xlsx.xlsx'
xml_files_folder = r'\\geneetika\Illumina_tsütokiip'
checked_files_log = "checked_xml_files.log"
CNV_fail_asukoht = r'\\srvlaste\Yhendlabor\GE_Illumina kiip\CNV_tabel_qSNP_GS_k6ik.bed'
# Testkeskkond
koondtabeli_asukoht = "GDA_38_koond.xlsx"
cnv_failide_kaust = r'C:\Users\Loom\Desktop\Testimiseks_2023_GDA'
# GS_faili_asukoht = r'\\srvlaste\Yhendlabor\GE_Illumina kiip\Kadi\GS_Kadi.xlsx'



# Siin algab qSNP failide läbivaatus

try:
    with open(checked_files_list,"r") as f:
        checked_files_list = [rida.rstrip('\r\n') for rida in list(f)]
except:
    checked_files_list = []

"""
def cnv_importija(asukoht):
    # print(asukoht)
    df = pd.read_csv(asukoht, sep='\t', lineterminator='\n')
    # Edasise sorteerimise tarbeks sundida tulp numbriks, "apply" ei sobi, sest algandmetes sees "date" kujul vigu
    df[['Max. Log BF','Chromosome']] = df[['Max. Log BF','Chromosome']].apply(pd.to_numeric, errors='coerce')
    # Duplikaatide kergemaks eristamiseks proovide ID-d ka sees
    df = df.loc[:,
        ['Sample Name','Chromosome',
        'Start Position (bp)','End Position (bp)',
        'Length (bp)','Copy Number','Max. Log BF',
        'Start Probe ID','End Probe ID']]
    usutavuse_filter = (df['Max. Log BF'] >= 10)
    #  ? return
    return(df.loc[usutavuse_filter])
"""

def process_xml_files(xml_files_folder, checked_files_list, checked_files_log):
    try:
        xml_file_locations = []
        for root, dirs, files in os.walk(xml_files_folder, topdown=False):
            if os.path.basename(root).startswith("GDA") and "Bookmark Analyses" in dirs:
                bookmark_folder = os.path.join(root, "Bookmark Analyses")
                for name in os.listdir(bookmark_folder):
                    if name.lower().endswith(".xml") and name.lower() not in checked_files_list:
                        xml_file_path = os.path.join(bookmark_folder, name)
                        xml_file_locations.append(xml_file_path)
                        with open(checked_files_log, 'a') as log_file:
                            log_file.write('\n' + xml_file_path.lower())
        # Update the list of checked files
        checked_files_list.extend([file.lower() for file in xml_file_locations])
    except Exception as e:
        print(f"An error occurred while processing XML files: {e}")

# Example usage:
cnv_files_folder = "/path/to/GDA_root_folder"
checked_files_list = []  # List to keep track of checked XML files
checked_files_log = "checked_xml_files.log"  # Log file to store checked XML file paths

process_xml_files(cnv_files_folder, checked_files_list, checked_files_log)

print("xml processed")




try:
    cnv_failide_asukohad = []
    for root, dirs, files in os.walk(cnv_failide_kaust, topdown=False):
        for name in files:
            if name.lower().endswith(".cnv") and name.lower() not in failide_nimekiri:
                cnv_failide_asukohad.append(os.path.join(root, name))
                with open(CNV_vaadatud_failid, 'a') as f:
                    f.write('\n' + name.lower())
except:
    print("Tekkis probleem cnv-failide asukohtade läbivaatamisel")

try:
    cnv_koondtabel = []
    if len(cnv_failide_asukohad) > 0:
        for cnv_faili_asukoht in cnv_failide_asukohad:
            yhe_patsiendi_cnvd = cnv_importija(cnv_faili_asukoht)
            print(f"Õnnestus import asukohast {cnv_faili_asukoht}")
            cnv_koondtabel.append(yhe_patsiendi_cnvd)
        cnv_koondtabel = pd.concat(cnv_koondtabel, ignore_index=True)
    else:
        print("Uusi cnv-faile ei ole")
except:
    print("Tekkis tõrge cnv-failidest koondtabeli tegemisel")

try:
    cnv_koondtabel['Length (bp)'] = cnv_koondtabel['Length (bp)'].apply(pd.to_numeric, errors='coerce')
    alam66duliste_filter = cnv_koondtabel['Length (bp)'] <= 30000
    cnv_koondtabel.drop(index=cnv_koondtabel[alam66duliste_filter].index)
    print("qSNP andmed imporditud")
except:
    print("Tekkis probleem qSNP cnv_koondtabeli pikkuse lahtri läbivaatamisel")



# Siin algab GenomeStudio andmete import

try:
    with open(GS_vaadatud_read,"r") as f:
        ridade_nimekiri = [rida.rstrip('\r\n') for rida in list(f)]
except:
    ridade_nimekiri = []
    print("GS vaadatud ridade nimekiri on tühi")
    
try:
    df_GS = pd.read_excel("GS_Kadi.xlsx")
    # Proovikoodi tagant on vaja ära saada [nr]
    df_GS['Sample ID'] = df_GS['Sample ID'].str[:9]
    vaadatud_ridade_filter = ~df_GS['Sample ID'].isin(ridade_nimekiri)
    df_GS = df_GS[vaadatud_ridade_filter]
except:
    print("Tekkis probleem GS_Kadi.xlsx faili lugemisel ja filtreerimisel")

try:
    # GS andmed kohandatud nii, et see ühilduks qSNP tabeliga
    df_GS = df_GS.loc[:,['Sample ID','Chr',
                         'Start','End',
                         'Size [bases]','Value','Comment']]
    df_GS = df_GS.rename(columns={'Sample ID':'Sample Name','Chr':'Chromosome',
                                  'Start':'Start Position (bp)','End':'End Position (bp)',
                                  'Size [bases]':'Length (bp)','Value':'Copy Number','Comment':'Max. Log BF'})
    # Sorteerimiseks on GS tabelis vaja X, Y ja XY (!) ära kaotada
    df_GS['Chromosome'] = df_GS['Chromosome'].replace({"XY":23})
    df_GS['Chromosome'] = df_GS['Chromosome'].replace({"X":23,"Y":24})
except:
    print("Tekkis probleem df_GS tulpade ümbernimetamisel ja ümberkodeerimisel")

# Genome Studio annab tohutult artefakte sugukromosoomides. Need tuleb maha rookida. Nende seas tõesed pikad CNV-d on qSNPs ka.
# Tüüpilised Genome Studio artefaktid on X-kromosoomis algusega 2786568 (peale PAR-regiooni) ja pikkusega 1,9-55 Mb; 
# ... Y-s levinud alguspunktid 2787139, 11500984 ja 18666909 ning pikkus Y-s 3,7-6 Mb
try:
    artefaktide_filter = ((df_GS['Start Position (bp)'] > 2786000) & (df_GS['Chromosome'] >= 23) & (df_GS['Length (bp)'] > 3700000))
    df_GS = df_GS[-artefaktide_filter]
except:
    print("Tekkis probleem df_GS sugukromosoomide artefaktide väljarookimisel")

try:
    # duplikaadid maha AINULT abifailis
    df_GS_IDd = df_GS['Sample Name'].drop_duplicates()
    with open(GS_vaadatud_read, 'a') as f:
        for ID in df_GS_IDd:
            f.write('\n' + ID)
    print('GS andmed imporditud')
except:
    print("Tekkis probleem GS duplikaatide eemaldamisel ja GS abifaili kirjutamisel")

# Siin lõppeb GenomeStudio andmete import


# Automaattuvastatud CNV-d exceli faili, kust loetakse hiljem sisse cnv-d juba läbi vaadatud qsnp failidest
try:
    vana_koondtabel = pd.read_excel("qSNP_k6ik_CNVd.xlsx")
    t2iendatud_koondtabel = pd.concat([vana_koondtabel, cnv_koondtabel, df_GS], ignore_index=True)
    praegune_aeg = datetime.now().strftime("%Y-%m-%d_%H-%M")
    os.rename("qSNP_k6ik_CNVd.xlsx",f"qSNP_k6ik_CNVd_backup_{praegune_aeg}.xlsx")
    df = t2iendatud_koondtabel.sort_values(by=['Chromosome','Start Position (bp)','End Position (bp)'])
    df.to_excel("qSNP_k6ik_CNVd.xlsx",index=False)
    print('Kõigi automaattuvastatud CNV-de Exceli fail uuendatud')
except:
    print("Tekkis probleem automaattuvastatud CNV-de Exceli faili uuendamisel")



# Siit algab BED-faili koostamine

try:
    os.rename("CNV_tabel_qSNP_GS_k6ik.bed",f"CNV_tabel_qSNP_GS_k6ik_backup_{praegune_aeg}.bed")
    # Alguses teha tühi fail vajaliku päisega. Sellele pannakse lõpuks DataFrame otsa
    with open("CNV_tabel_qSNP_GS_k6ik.bed", 'a') as f:
        f.write('track db="hg38" name="GenK k6ik CNVd" description="Geneetikakeskuse qSNP ja Genome Studio CNV-d. GRCh38" itemRgb="On"' + '\n')
    # Duplikaatide esile tõstmiseks saaks kasutada alguse ja lõpu proovide ID tulpi
    df = df.loc[:,['Chromosome','Start Position (bp)','End Position (bp)',
                   'Sample Name','Length (bp)','Copy Number','Start Position (bp)','End Position (bp)']]
except:
    print("Tekkis probleem vana qSNP BED-faili ümbernimetamisel ja uue põhja loomisel")

try:
    # Plaadi nr ja muu lisainfo tuleb GRCh38 koondtabelist serveris
    df2 = pd.read_excel(koondtabeli_asukoht)
    df2 = df2.rename(columns={'SampleID eLabor':'Sample Name'})
except:
    print("Tekkis probleem BED-faili tooriku eeltöös: GRCh38 leidude koondtabeli impordil")


def sobimatute_CNVde_eemaldaja (vastusekuupaev,aberratsioon):
    # Eemaldada kehva DNA kvaliteediga artefakte täis proovid ("korrata","tühistatud","analüüsimatu")
    if "kor" in str(vastusekuupaev) or "tüh" in str(vastusekuupaev) or "analüüsima" in str(vastusekuupaev):
        return None
    # Eemaldada triploidiad ja trisoomiad, sugukromosoomide anomaaliad las jäävad - autosoomide CNVd las olla näha. Hulgileiud jäävad.
    if ")x" in str(aberratsioon) or "trip" in str(aberratsioon) or "kontam" in str(aberratsioon):
        return None
    else:
        return "korras"
    
    
try:
    df2['Sobib'] = df2.apply(lambda x: sobimatute_CNVde_eemaldaja(x.Vastus, x.Aberratsioon), axis=1)
except:
    print("Tekkis probleem sobimatute CNV-de eemaldamisel BED-faili toorikust")
    
try:
    df = pd.merge(df,df2[['Sample Name','Plaat','Patsient','Sobib']],how='left',on='Sample Name')
    df = df.dropna(subset=['Sobib'])
    df[['Chromosome','Plaat','Copy Number']] = df[['Chromosome','Plaat','Copy Number']].fillna(0.0).astype(int)
except:
    print("Tekkis probleem BED-tooriku ja CRCh38 koondtabeli ühendamisel")


# Tavapatsientidest eristada hematoloogilised ja emad-isad
def patsiendi_tyybiti_sorteerija(patsient):
    if patsient[:3] == 'loo':
        return 'loode/'
    if patsient[:3] == 'hem':
        return 'hematol./'
    if patsient[:3] == 'ema' or patsient[:3] == 'isa':
        return 'vanem/'    
    else:
        return ''

try:
    df['Patsient'] = df['Patsient'].astype(str).apply(patsiendi_tyybiti_sorteerija)
    df['Sample Name'] = df['Patsient'] + df['Plaat'].astype(str) + '/' + df['Sample Name']
    df = df.drop(columns=['Plaat','Patsient','Sobib'])
except:
    print("Tekkis probleem BED-toorikus patsientide tüübiti sorteerimisel")


# Koopia-arvu lahtris tohivad olla ainult numbrid
def koopiaarv_v2rviks(copy_number):
    if copy_number >= 3:
        return '0,0,255'
    elif copy_number <= 1:
        return '255,0,0'
    elif copy_number == 2:
        return '150,50,255'
    else:
        pass
        

try:
    df['BED_v2rv'] = df['Copy Number'].apply(koopiaarv_v2rviks)
    # nüüdseks ebavajalikest tulpadest saavad BED-faili kohustuslikud tühjad tulbad
    df.loc[:,'Length (bp)'] = 0
    df.loc[:,'Copy Number'] = "."
    # Asendus 23->X tuleb teha lõpus, et sorteerimisvõimekust säilitada
    df['Chromosome'] = df['Chromosome'].replace({23:"X",24:"Y"})
    # kogu tulp vaja stringiks teha, muidu ei saa "chr" juurde kirjutada
    df['Chromosome'] = df['Chromosome'].astype(str)
    df['Chromosome'] = "chr" + df['Chromosome']
except:
    print("Tekkis probleem BED-tooriku värvi-, kromosoomi- ja abilahtrite koostamiselt")
    
# try:
#    df.to_csv(CNV_fail_asukoht, header=0, index=None, sep='\t', mode='a')
#except:
#    print("Tekkis probleem teises kaustas oleva qSNP/GS BED-faili loomisega")
    
try:
    df.to_csv(r'CNV_tabel_qSNP_GS_k6ik.bed', header=0, index=None, sep='\t', mode='a')
    print('Valmis!')
except:
    print("Tekkis probleem uue BED-faili kirjutamisel")