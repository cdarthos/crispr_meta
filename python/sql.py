import pandas as pd

import function
import requests
import flask
from Bio import Entrez
import xml.etree.ElementTree as ET
import os.path as p
import xmltodict
import sqlite3
import pandas as pb
import csv
import sys



class bdd_sqlite:
    def __init__(self, db_name=':memory:'):
        self.db_name = db_name
        self.conn = sqlite3.connect(self.db_name)
        self.cursor = self.conn.cursor()
        #self.create_db_ncbi()
        self.create_db_3651ref()
        self.extract_all_3651()

    # cree la table ncbi
    def create_db_ncbi(self, sql="create_table_ncbi.sql"):
        with open(sql) as request:
            sql_request = request.read()
            self.cursor.executescript(sql_request)

    # cree la table de corrspondance entre Isolate accession number & ncbi
    def create_db_ebi2ncbi(self, sql="create_table_ebi2ncbi.sql"):
        with open(sql) as request:
            sql_request = request.read()
            self.cursor.executescript(sql_request)

    # cree la table de 3651ref
    def create_db_3651ref(self):
        # self.cursor.execute('''CREATE TABLE ref_3651 (
        #                                                         Isolate_accession_number text,
        #                                                         Isoniazid text,
        #                                                         Rifampicin text,
        #                                                         Ethambutol text,
        #                                                         Pyrazinamide text,
        #                                                         Streptomycin text,
        #                                                         Ciprofloxacin text,
        #                                                         Moxifloxacin text,
        #                                                         Ofloxacin text,
        #                                                         Amikacin text,
        #                                                         Capreomycin text,
        #                                                         Kanamycin text)''')
        imp_3651_ref = pd.read_csv('SRAlist/3651ref.csv',sep=';', encoding='latin-1')
        imp_3651_ref.to_sql('ref_3651',self.conn,if_exists='append', index = False)

        #df = pd.read_sql('SELECT * FROM ref_3651 WHERE ref_3651.Isolate_accession_number="ERS456778"', self.conn)
        # Verify that result of SQL query is stored in the dataframe
        #print(df.head())

    def extract_all_3651(self):
        with open("SRAlist/ebi_2_run_acc_file.csv", 'r', encoding='utf-8') as input:
            acc_key_value = input.readline().replace("\n", "")



            acc_key_value = input.readline().replace("\n", "")
            acc_key_value = acc_key_value.split(";")
            dict_ncbi = function.find_sra_element(acc_key_value[1])
            sql_req = 'SELECT * FROM ref_3651 WHERE ref_3651.Isolate_accession_number="' + acc_key_value[
                0] + '" LIMIT 1'
            dict_ebi = pd.read_sql(sql_req, self.conn).to_dict(orient='record')
            dict_final = {**dict_ncbi, **dict_ebi[0]}
            header = dict_final.keys()
            with open("SRAlist/ebiANDncbi.csv", 'w', encoding='utf-8', newline='') as out:

                writer = csv.DictWriter(out, delimiter='\t', fieldnames=header)
                writer.writeheader()

                while out != "":
                    acc_key_value = input.readline().replace("\n", "")
                    acc_key_value = acc_key_value.split(";")
                    dict_ncbi = function.find_sra_element(acc_key_value[1])
                    sql_req = 'SELECT * FROM ref_3651 WHERE ref_3651.Isolate_accession_number="' + acc_key_value[0] +'"'
                    dict_ebi = pd.read_sql(sql_req, self.conn).to_dict(orient='records')
                    try:
                        #print(dict_ebi[0])
                        dict_final = {**dict_ncbi, **dict_ebi[0]}
                        writer.writerow(dict_final)
                    except:
                        print(dict_ebi)
                        print(acc_key_value)
                        print("Unexpected error:", sys.exc_info())
                        print(acc_key_value[0] + " / " + acc_key_value[1])






    # Ajoute dans la table ncbi les informations issue de NCBI pour chaque SRA présent dans le fichier
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
                self.insert_table_ncbi(
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

    def insert_table_ncbi(self, task):
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


if __name__ == "__main__":
    bd = bdd_sqlite()
