import csv

import flask
from Bio import Entrez
import xml.etree.ElementTree as ET
import os.path as p
import xmltodict, json

app = flask.Flask(__name__)

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0')


@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"


@app.route("/<Run_accession>")
def json_sra(Run_accession):
    fetch_sra_fullXML(Run_accession)
    with open("complete_xml/" + Run_accession + ".xml", 'r') as xml:
        o = xmltodict.parse(xml.read())
        #return o.get("BioSample")
        return o


def fetch_sra_fullXML(id):
    outname = id + ".xml"
    if p.exists("complete_xml/" + outname) and p.getsize("complete_xml/" + outname) != 0:
        return outname
    Entrez.email = "clement.lecarpentier02@edu.univ-fcomte.fr"
    Entrez.max_tries = 5
    Entrez.sleep_between_tries = 30
    handle = Entrez.efetch(db="sra", id=id, rettype='full', retmode="xml")
    # rec1 = Entrez.read(handle)
    # print(handle)
    record = handle.read().decode('UTF-8').replace('\n', '')
    handle.close()
    with open("complete_xml/" + outname, 'w', encoding="utf-8") as out:
        out.write(record)
    return id


def print_sra(Run_accession):
    fetch_sra_fullXML(Run_accession)
    with open("complete_xml/" + Run_accession + ".xml", 'r') as xml:
        tree = ET.parse(xml)
        root = tree.getroot()
        # print('Item #2 attribute:')
        # print(root[0][1].attrib)
        for user in root.iter():
            print(user.tag, user.items(), user.text)


def fetch_BioProject_fullXML(id):
    outname = id + ".xml"
    if (p.exists("BioProject/" + outname)):
        return id
    Entrez.email = "clement.lecarpentier02@edu.univ-fcomte.fr"
    Entrez.max_tries = 5
    Entrez.sleep_between_tries = 30
    handle = Entrez.efetch(db="BioProject", id=id, rettype='full', retmode="xml")
    # rec1 = Entrez.read(handle)
    # print(handle)
    record = handle.read().decode('UTF-8').replace('\n', '')
    handle.close()
    with open("BioProject/" + outname, 'w', encoding="utf-8") as out:
        out.write(record)
    return id


def print_BioProject(path):
    with open("BioProject/" + path + ".xml", 'r') as xml:
        tree = ET.parse(xml)
        root = tree.getroot()
        # print('Item #2 attribute:')
        # print(root[0][1].attrib)

        for user in root.iter():
            print(user.tag, user.items(), user.text)


def file_sra_req(file):
    with open("SRAlist/" + file, 'r') as SRAfile:
        max = 400
        i = 0
        SRAlist = SRAfile.readline()
        while SRAlist != "":
            # print(SRAlist)
            out = fetch_sra_fullXML(SRAlist.strip())
            print(str(i) + "/87122")
            i += 1
            # print(out)
            SRAlist = SRAfile.readline()

@app.route("/extract/<Run_accession>")
def find_sra_element(Run_accession):
    fetch_sra_fullXML(Run_accession)
    keys_list = ['Run_accession', 'datePubli', 'sraID', 'BioSample','collection_date', 'geo_loc_name',
                 'lat_lon','list_extra','BioProject', 'title', 'abstract']
    resultat_dict = dict.fromkeys(keys_list)
    resultat_dict['Run_accession'] = Run_accession
    with open("complete_xml/" + Run_accession + ".xml", 'r', encoding="utf-8", errors='ignore') as xml:
        tree = ET.parse(xml)
        root = tree.getroot()
        for n in root.findall(".//STUDY/DESCRIPTOR/STUDY_TITLE"):
            resultat_dict['title'] = n.text
            # print(".//STUDY/DESCRIPTOR/STUDY_TITLE",n.text)
        for n in root.findall(".//STUDY/DESCRIPTOR/STUDY_ABSTRACT"):
            resultat_dict['abstract'] = n.text
            # print(".//STUDY/DESCRIPTOR/STUDY_ABSTRACT",n.text)
        # for n in root.findall(".//SUBMISSION"):

        # print(".//SUBMISSION",n.items())
        for n in root.findall(".//SUBMISSION"):
            n = n.get('accession')
            resultat_dict['sraID'] = n
            # print(".//SUBMISSION:accession",n)

        for n in root.findall(".//SCIENTIFIC_NAME"):
            resultat_dict['specie'] = n.text
            # print(".//SAMPLE/SAMPLE_NAME/SCIENTIFIC_NAME",n.text)
        for n in root.findall(".//IDENTIFIERS/EXTERNAL_ID[@namespace='BioSample']"):
            if(n.text != ""):
                resultat_dict['BioSample'] = n.text
            # print("biosampleLABEL = ",n.text)
        for n in root.findall(".//EXTERNAL_ID[@namespace='BioProject']"):
            if (n.text != ""):
                resultat_dict['BioProject'] = n.text
            # print("biosampleLABEL = ",n.text)
        list_extra = []
        for n in root.findall(".//SAMPLE_ATTRIBUTE"):
            extra = {}
            for o in n.iter():
                if (o.tag == "SAMPLE_ATTRIBUTE"):
                    continue
                elif (o.tag == "TAG"):
                    tag = o.text
                    continue
                else:
                    extra[tag] = o.text
                    # print(tag, o.text)
            list_extra.append(extra)
        resultat_dict['list_extra'] = list_extra
        for extra in list_extra:
            if([*extra] == ['geo_loc_name']):
                #print([*extra.values()][0])
                resultat_dict['geo_loc_name'] = [*extra.values()][0]
            if([*extra] == ['lat_lon']):
                #print([*extra.values()][0])
                resultat_dict['lat_lon'] = [*extra.values()][0]
            if([*extra] == ['collection_date']):
                #print([*extra.values()][0])
                resultat_dict['collection_date'] = [*extra.values()][0]

        for n in root.findall(".//RUN_SET/RUN"):
            resultat_dict['datePubli'] = n.get("published")
            # print(".//RUN_SET/RUN:published",datePubli)
        # for n in root.iter():
        # print(n.tag, n.attrib, n.text)
    return resultat_dict



def extract_in_csv(file):
    with open("SRAlist/" + file, 'r') as SRAfile:
        i = 0
        # header = ['Run_accession', 'datePubli', 'pubmed_id', 'pubmedLABEL', 'sraID', 'biosampleLABEL', 'country', 'list_extra', 'title', 'abstract']

        SRAlist = SRAfile.readline()

        with open("test.csv", 'w', newline='', encoding='utf-8') as out:
            resultat_dict = find_sra_element(SRAlist.strip())
            header = resultat_dict.keys()
            # header = ['Run_accession', 'datePubli', 'pubmed_id', 'pubmedLABEL', 'sraID', 'biosampleLABEL', 'country', 'list_extra', 'title', 'abstract']
            # print(header)

            writer = csv.DictWriter(out, delimiter='\t', fieldnames=header)
            writer.writeheader()
            # header = "Run_accession\tdatePubli\tpubmed_id\tpubmedLABEL\tsraID\tbiosampleLABEL\tcountry\tlist_extra\ttitle\tabstract"
            # out.write(header + "\n")
            while SRAlist != "":
                out = fetch_sra_fullXML(SRAlist.strip())
                print(SRAlist + " => " + str(i) + "/87122")
                i += 1

                resultat_dict = find_sra_element(SRAlist.strip())
                print(resultat_dict)
                writer.writerow(resultat_dict)

                # line = [Run_accession, datePubli, pubmed_id, pubmedLABEL, sraID, biosampleLABEL, country, list_extra, title, abstract]
                # line = [title, abstract]

                # writer.writerow(line)

                SRAlist = SRAfile.readline()

@app.route("/update")
def update():
    #extract_in_csv("SRAlistTEST.txt")
    file_sra_req("SRAlistTEST.txt")

@app.route("/update_csv")
def update_csv():
    extract_in_csv("SRAlistTEST.txt")
