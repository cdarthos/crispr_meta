import function

import requests

import flask
from Bio import Entrez
import xml.etree.ElementTree as ET
import os.path as p
import xmltodict
import sqlite3

from flask import jsonify

app = flask.Flask(__name__)


def ebi_2_run_acc_file():
    with open("SRAlist/ebi", 'r') as SRAfile:
        with open("SRAlist/ebi_2_run_acc_file.csv", 'w', encoding='utf-8') as out:
            out.write("id_ebi;Run_accession\n")
        SRAlist = SRAfile.readline()
        id_ebi = SRAlist.split(",")
        for id in id_ebi:
            run_acc_list = ebi_2_run_acc(id.replace("\n",""))
            run_acc_list = run_acc_list.split("\n")
            for run_acc in run_acc_list:
                if (run_acc != ""):
                    with open("SRAlist/ebi_2_run_acc_file.csv", 'a', encoding='utf-8') as out:
                        out.write(id.replace("\n","") + ";" + run_acc + "\n")

        while SRAlist != "":
            SRAlist = SRAfile.readline()
            id_ebi = SRAlist.split(",")
            for id in id_ebi:

                run_acc_list = ebi_2_run_acc(id.replace("\n", ""))
                run_acc_list = run_acc_list.split("\n")
                for run_acc in run_acc_list:
                    if run_acc:
                        print(id.replace("\n", "") + ";" + run_acc)
                        with open("SRAlist/ebi_2_run_acc_file.csv", 'a', encoding='utf-8') as out:
                            out.write(id.replace("\n", "") + ";" + run_acc + "\n")

    return "OK"

def ebi_2_run_acc(ebi_id):
    req = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + ebi_id + "&result=read_run&fields=run_accession"
    r = requests.get(req)
    run_acc_temp = r.text
    run_acc = run_acc_temp.replace("run_accession\n", "")
    return run_acc


if __name__ == "__main__":
    #app.run(host="0.0.0.0", debug=True)
    out = ebi_2_run_acc_file()



@app.route("/")
def hello_world():
    return "<p>Hello, World!</p>"


@app.route("/<Run_accession>")
def json_sra(Run_accession):
    function.fetch_sra_fullXML(Run_accession)
    with open("complete_xml/" + Run_accession + ".xml", 'r') as xml:
        o = xmltodict.parse(xml.read())
        return o


@app.route("/update")
def update():
    # extract_in_csv("SRAlistTEST.txt")
    function.file_sra_req("SRAlistTEST.txt")


@app.route("/update_csv")
def update_csv():
    function.extract_in_csv("SRAlistTEST.txt")
    return ("OK")


@app.route("/extract/<Run_accession>")
def find_sra_element_f(Run_accession):
    out = function.find_sra_element(Run_accession)
    return out


@app.route("/bd")
def bdd_flask():
    bd = bdd_sqlite()
    out = bd.import_sra()
    return str(out)


@app.route("/ebi/<ebi_id>")
def ebi_2_run_acc(ebi_id):
    req = "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=" + ebi_id + "&result=read_run&fields=run_accession"
    r = requests.get(req)
    run_acc_temp = r.text
    run_acc = run_acc_temp.replace("run_accession\n", "")
    return run_acc


@app.route('/ebi/test')
def ebi_2_run_acc_file():
    with open("SRAlist/ebi", 'r') as SRAfile:
        with open("SRAlist/ebi_2_run_acc_file.csv", 'w', encoding='utf-8') as out:
            out.write("id_ebi;Run_accession\n")
            SRAlist = SRAfile.readline()
            id_ebi = SRAlist.split(",")
            for id in id_ebi:
                run_acc_list = ebi_2_run_acc(id).split("\n")
                for run_acc in run_acc_list:
                    if (run_acc != ""):
                        out.write(id + ";" + run_acc + "\n")
                        print(id + ";" + run_acc + "\n")
            while SRAlist != "":
                SRAlist = SRAfile.readline()
                id_ebi = SRAlist.split(",")
                for id in id_ebi:
                    run_acc_list = ebi_2_run_acc(id).split("\n")
                    for run_acc in run_acc_list:
                        if (run_acc != ""):
                            out.write(id + ";" + run_acc + "\n")
                            print(id + ";" + run_acc + "\n")
    return "OK"

class bdd_sqlite:
    def __init__(self, db_name=':memory:'):
        self.db_name = db_name
        self.conn = sqlite3.connect(self.db_name)
        self.cursor = self.conn.cursor()
        self.create_db()

    def import_sra(self, file="test1.txt"):
        i = 0
        with open("SRAlist/" + file, 'r') as SRAfile:
            # SRAlist = SRAfile.readline()
            # dict = function.find_sra_element(SRAlist.strip())
            SRAlist = "OK"
            while SRAlist != "":
                SRAlist = SRAfile.readline()
                dict = function.find_sra_element(SRAlist.strip())
                i += 1
                print("SRA numero : ", i)
                self.insert_table(
                    [dict["Run_accession"],
                     dict["BioSample"],
                     dict["BioProject"],
                     dict["collection_date"],
                     dict["geo_loc_name"],
                     dict["isolate"],
                     dict["lat_lon"],
                     dict["sraID"],
                     dict["strain"],
                     dict["title"]]
                )
            return i

    def create_db_ncbi(self, sql="create_table_ncbi.sql"):
        with open(sql) as request:
            sql_request = request.read()
            self.cursor.executescript(sql_request)

    def show_table(self):
        sql = '''
                SELECT *
                FROM ncbi 
                ;
        '''
        try:
            self.cursor.execute(sql)
        except sqlite3.Error as e:
            print(e)
        list_dict = []
        for line in self.cursor.fetchall():
            dict = {}
            dict["Run_accession"] = line[0]
            dict["BioSample"] = line[1]
            dict["BioProject"] = line[2]
            dict["collection_date"] = line[3]
            dict["geo_loc_name"] = line[4]
            dict["isolate"] = line[5]
            dict["lat_lon"] = line[6]
            dict["sraID"] = line[7]
            dict["strain"] = line[8]
            dict["title"] = line[9]
            list_dict.append(dict)
        return list_dict

    def insert_table(self, task):
        sql = '''INSERT INTO ncbi (
                        Run_accession ,
                        BioSample ,
                        BioProject ,
                        collection_date ,
                        geo_loc_name ,
                        isolate ,
                        lat_lon ,
                        sraID ,
                        strain,
                        title
                        )
                    VALUES (?,?,?,?,?,?,?,?,?,?)        
        '''
        self.cursor.execute(sql, task)
        self.conn.commit()

    def close(self):
        self.conn.close()
