# Tool to combine all xml files from directory to single BED file.
# Input: Directory with xml files, patient data from Excel file, text files to keep track of sample numbers
# Output: xlsx-file and BED-fail of all automatically detected CNV-s.
# 29.02.2024
# v4.2 (GDA int format fix, without .0)
# Kadi Jairus


import os
import pandas as pd
from datetime import datetime
import xml.etree.ElementTree as ET

print("Programm alustas tööd, palun oota...")

# Files in server
GDA_patients_table = r'\\srvlaste\Yhendlabor\GE_Illumina kiip\GDA_38_koond.xlsx'
xml_files_folder = r'\\geneetika\Illumina_tsütokiip'
checked_files_log = "checked_xml_files.log"
CNV_file_location = r'\\srvlaste\Yhendlabor\GE_Illumina kiip\CNV_tabel_qSNP_GS_k6ik.bed'
# Files in test environment
# GDA_patients_table = "GDA_38_koond.xlsx"
# xml_files_folder = r'C:\Users\Admin\Desktop\Testimiseks_2023_GDA'
# xml_files_folder = r'C:\Users\Loom\Desktop\Testimiseks_2023_GDA'
CNV_file_location = "GDA_CNV_tabel_GS_k6ik.bed"


try:
    with open(checked_files_log,"r") as f:
        checked_files_list = [rida.rstrip('\r\n') for rida in list(f)]
        print(checked_files_list)
except:
    print("Läbivaadatud failide nimekiri on tühi.")
    checked_files_list = []

def process_xml_files(xml_files_folder, checked_files_list, checked_files_log):
    try:
        # print(xml_files_folder)
        xml_file_locations = []
        for root, dirs, files in os.walk(xml_files_folder, topdown=False):
            if any(os.path.basename(d).startswith("GDA") for d in root.split(os.path.sep)) and "Bookmark Analyses" in dirs:
                bookmark_folder = os.path.join(root, "Bookmark Analyses")
                for name in os.listdir(bookmark_folder):
                    if name.lower().endswith(".xml") and name.lower() not in checked_files_list:
                        xml_file_path = os.path.join(bookmark_folder, name)
                        xml_file_locations.append(xml_file_path)
                        with open(checked_files_log, 'a') as log_file:
                            log_file.write('\n' + xml_file_path.lower())
        # Update the list of checked files
        checked_files_list.extend([file.lower() for file in xml_file_locations])
        return xml_file_locations
    except Exception as e:
        print(f"An error occurred while processing XML files: {e}")


xml_file_locations = process_xml_files(xml_files_folder, checked_files_list, checked_files_log)

# print("xml processed")
# print(xml_file_locations)


def import_cnv_data_from_xml(xml_file_location):
    """Function to import CNV data from an XML file into a Pandas DataFrame."""
    tree = ET.parse(xml_file_location)
    root = tree.getroot()

    data = []
    for bookmark_elem in root.findall('.//bookmark'):
        sample_id = bookmark_elem.find('sample_id').text[:9]
        chr_num = bookmark_elem.find('chr_num').text
        # Replace 'XY' with 'X' in PAR region.
        if chr_num == 'XY':
            chr_num = 'X'
        base_start_pos = int(bookmark_elem.find('base_start_pos').text)
        base_end_pos = int(bookmark_elem.find('base_end_pos').text)
        value = int(bookmark_elem.find('value').text)

        size = base_end_pos - base_start_pos
        if size >= 30000:
            data.append([sample_id, chr_num, base_start_pos, base_end_pos, size, value])

    # columns = ['sample_id', 'chr_num', 'base_start_pos', 'base_end_pos', 'size', 'value']
    columns = ['Sample Name', 'Chromosome', 'Start Position (bp)', 'End Position (bp)', 'Length (bp)', 'Copy Number']
    df = pd.DataFrame(data, columns=columns)
    return df


def extract_data_from_xml_file_locations(xml_file_locations):
    """Using the locations of newest xml-file(s), get XML to DataFrame."""
    cnvs_from_xml = []

    try:
        for xml_file_location in xml_file_locations:
            cnvs_of_single_patient = import_cnv_data_from_xml(xml_file_location)
            print("Edukas import")
            print(f"Successfully imported data from {xml_file_location}")
            cnvs_from_xml.append(cnvs_of_single_patient)
            print(cnvs_from_xml)
        cnvs_from_xml = pd.concat(cnvs_from_xml, ignore_index=True)
        return cnvs_from_xml
    except Exception as e:
        print(f"An error occurred while processing CNV files: {e}")
        return None


# Define the custom sorting order for chr_num
chr_num_order = [
    '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12',
    '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X', 'Y'
]


def append_new_to_old_cnv_data(df_new_cnvs):
    try:
        # Read existing CNVs from Excel file
        df_old_cnvs = pd.read_excel("Genome_Studio_k6ik_CNVd.xlsx")
        df_cnvs_appended = pd.concat([df_old_cnvs, df_new_cnvs], ignore_index=True)
        # print(df_cnvs_appended)
        current_time = datetime.now().strftime("%Y-%m-%d_%H-%M")
        os.rename("Genome_Studio_k6ik_CNVd.xlsx", f"Genome_Studio_k6ik_CNVd_backup_{current_time}.xlsx")
        # Convert chr_num to categorical data with custom ordering.
        df_cnvs_appended['Chromosome'] = pd.Categorical(df_cnvs_appended['Chromosome'], categories=chr_num_order, ordered=True)
        df_sorted = df_cnvs_appended.sort_values(by=['Chromosome', 'Start Position (bp)', 'End Position (bp)'])
        df_sorted.to_excel("Genome_Studio_k6ik_CNVd.xlsx", index=False)
        print('Kõigi automaattuvastatud CNV-de Exceli fail uuendatud')
        return df_sorted
    except Exception as e:
        print(f"Tekkis probleem automaattuvastatud CNV-de Exceli faili uuendamisel: {e}")

# Example usage:
df_new_cnvs = extract_data_from_xml_file_locations(xml_file_locations)

# Check if the result is not None before further processing
if df_new_cnvs is not None:
    # Continue with further processing or analysis
    print("Uued CNV-d on:")
    print(df_new_cnvs.head(100))
    df_all_cnvs_sorted = append_new_to_old_cnv_data(df_new_cnvs)

def remove_duplicate_cnvs(df):
    """Keep only unique CNV-s of the same patient."""
    # Find duplicates based on specified columns
    duplicate_rows = df[df.duplicated(subset=['Sample Name', 'Chromosome', 'Start Position (bp)', 'End Position (bp)'], keep=False)]

    # Count number of duplicate rows
    num_duplicates = len(duplicate_rows)

    # Remove duplicates and keep the first occurrence
    df_unique = df.drop_duplicates(subset=['Sample Name', 'Chromosome', 'Start Position (bp)', 'End Position (bp)'], keep='first')
    print(f"Number of removed duplicates: {num_duplicates}")
    return df_unique


def sobimatute_CNVde_eemaldaja (vastusekuupaev,aberratsioon):
    # Eemaldada kehva DNA kvaliteediga artefakte täis proovid ("korrata","tühistatud","analüüsimatu")
    if "kor" in str(vastusekuupaev) or "tüh" in str(vastusekuupaev) or "analüüsima" in str(vastusekuupaev):
        return None
    # Eemaldada triploidiad ja trisoomiad, sugukromosoomide anomaaliad las jäävad - autosoomide CNVd las olla näha. Hulgileiud jäävad.
    if ")x" in str(aberratsioon) or "trip" in str(aberratsioon) or "kontam" in str(aberratsioon):
        return None
    else:
        return "korras"


def add_patient_data_to_cnv(df):
    """Add patient data and filter."""
    try:
        """Get additional patient data from Excel"""
        df2 = pd.read_excel(GDA_patients_table)
        # Convert chr_num to categorical data with custom ordering.
        df['Chromosome'] = pd.Categorical(df['Chromosome'], categories=chr_num_order, ordered=True)
        print(df2)
        df2.to_csv("Vaheseis_GDA_patsiendid.csv", header=1, index=None, sep='\t',mode='a')
        df2 = df2.rename(columns={'SampleID eLabor':'Sample Name'})
    except Exception as e:
        print(f"Tekkis probleem BED-faili tooriku eeltöös: GDA koondtabeli impordil: {e}")

    try:
        df2['Sobib'] = df2.apply(lambda x: sobimatute_CNVde_eemaldaja(x.Vastuse_kp, x.Aberratsioon), axis=1)
    except Exception as e:
        print(f"Tekkis probleem probleemsete CNV-de eemaldamisel: {e}")

    try:
        df = pd.merge(df,df2[['Sample Name','Plaat','Patsient','Sobib']],how='left',on='Sample Name')
        # df = df.dropna(subset=['Sobib'])
        df.dropna(subset=['Sobib', 'Chromosome', 'Plaat', 'Copy Number'], inplace=True)
        # df[['Chromosome','Plaat','Copy Number']] = df[['Chromosome','Plaat','Copy Number']].fillna(0.0).astype(int)
        return df
    except Exception as e:
        print(f"Tekkis probleem BED-tooriku ja GDA koondtabeli ühendamisel: {e}")


def shorten_patient_type(patient):
    """Add label if patient is a (supposedly normal) parent or cancer patient."""
    if patient[:3] == 'loo':
        return 'loode/'
    if patient[:7] == 'hematol' or patient[:7] == 'tahke k':
        return 'onko/'
    if patient[:3] == 'ema' or patient[:3] == 'isa':
        return 'vanem/'    
    else:
        return ''

def add_label_to_special_patient_types(df):
    """Add label if patient is a (supposedly normal) parent or cancer patient."""
    try:
        df['Patsient'] = df['Patsient'].astype(str).apply(shorten_patient_type)
        df['Plaat'] = df['Plaat'].apply(lambda x: str(int(x)) if isinstance(x, (int, float)) else x)
        df['Sample Name'] = df['Patsient'] + "GDA" + df['Plaat'] + '/' + df['Sample Name']
        df = df.drop(columns=['Plaat','Patsient','Sobib'])
        
        #df['Patsient'] = df['Patsient'].astype(str).apply(shorten_patient_type)
        #df['Sample Name'] = df['Patsient'] + "GDA" + df['Plaat'].astype(str).replace(".0","") + '/' + df['Sample Name']
        #df = df.drop(columns=['Plaat','Patsient','Sobib'])
        return df
    except Exception as e:
        print(f"Tekkis probleem BED-toorikus patsientide tüübiti sorteerimisel: {e}")
        return None


def make_bed_file_from_df(df):
    """Create BED-file from filtered all CNV-s dataframe. """
    try:
        praegune_aeg = datetime.now().strftime("%Y-%m-%d_%H-%M")
        os.rename("GDA_CNV_tabel_GS_k6ik.bed",f"GDA_CNV_tabel_GS_k6ik_backup_{praegune_aeg}.bed")
        # Alguses teha tühi fail vajaliku päisega. Sellele pannakse lõpuks DataFrame otsa
        with open("GDA_CNV_tabel_GS_k6ik.bed", 'a') as f:
            f.write('track db="hg38" name="GenK k6ik GDA CNVd" description="Geneetikakeskuse automaattuvastatud GDA Genome Studio CNV-d. GRCh38" itemRgb="On"' + '\n')
        # Duplikaatide esile tõstmiseks saaks kasutada alguse ja lõpu proovide ID tulpi
        df = df.loc[:,['Chromosome','Start Position (bp)','End Position (bp)',
                       'Sample Name','Length (bp)','Copy Number','Start Position (bp)','End Position (bp)']]
        return df
    except Exception as e:
        print(f"Tekkis probleem vana BED-faili ümbernimetamisel ja uue põhja loomisel {e}")
        return None

df = remove_duplicate_cnvs(df_all_cnvs_sorted)
df_with_patient_data = add_patient_data_to_cnv(df)
print("df_with_patient_data:")
print(df_with_patient_data.head(50))
# df_with_labeled_patient_data
df = add_label_to_special_patient_types(df_with_patient_data)
print("df_with_labeled_patient_data")
print(df.head(50))

df = make_bed_file_from_df(df)


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
    # df['Chromosome'] = df['Chromosome'].replace({23:"X",24:"Y"})
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
    df.to_csv(r'GDA_CNV_tabel_GS_k6ik.bed', header=0, index=None, sep='\t', mode='a')
    print('Valmis!')
except:
    print("Tekkis probleem uue BED-faili kirjutamisel")