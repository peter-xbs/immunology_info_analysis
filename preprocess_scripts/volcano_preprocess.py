# _*_ coding:utf-8 _*_

import pickle
import requests
import os


# url = 'https://www.ncbi.nlm.nih.gov/geo/geo2r/backend/geo2r.cgi?ctg_time=1539095023&job_key=JSID_01_592643_130.14.22.10_9000_geo2r'
#
# url = 'https://www.ncbi.nlm.nih.gov/geo/geo2r/backend/geo2r.cgi?ctg_time=1541686008&job_key=JSID_01_614458_130.14.22.10_9000_geo2r'
def get_data(url, home_dir, gse_id):
    res = requests.get(url).content.decode()
    output = os.path.join(home_dir, gse_id+'.diff')
    with open(output, 'w') as fo:
        fo.write(res)

#url = "https://www.ncbi.nlm.nih.gov/geo/geo2r/backend/geo2r.cgi?ctg_time=1541845165&job_key=JSID_01_766447_130.14.18.128_9000_geo2r"
#url = "https://www.ncbi.nlm.nih.gov/geo/geo2r/backend/geo2r.cgi?ctg_time=1541864801&job_key=JSID_01_726343_130.14.18.6_9000_geo2r"
#url = "https://www.ncbi.nlm.nih.gov/geo/geo2r/backend/geo2r.cgi?ctg_time=1541902709&job_key=JSID_01_726477_130.14.18.6_9000_geo2r"
#url = "https://www.ncbi.nlm.nih.gov/geo/geo2r/backend/geo2r.cgi?ctg_time=1541904801&job_key=JSID_01_615897_130.14.22.10_9000_geo2r"
url = "https://www.ncbi.nlm.nih.gov/geo/geo2r/backend/geo2r.cgi?ctg_time=1541905198&job_key=JSID_01_726490_130.14.18.6_9000_geo2r"
key1 = "glycolytic process" #糖酵解
key2 = "gluconeogenesis" # 糖异生
key3 = "lipid metabolic process" # 脂代谢
key4 = "cellular amino acid metabolic process" # 氨基酸代谢过程
key5 = "carbohydrate metabolic process" # 碳酸代谢过程
key6 = "glutamate metabolic process " # 谷氨酸代谢
key7 = "glucose metabolic process" # 葡萄糖代谢过程
key8 = "reactive oxygen species metabolic process" # ROS 代谢
key9 = "one-carbon metabolic process" # 单碳循环
key10 = "glutamine metabolic process" #
key11 = "NAD metabolic process"
key12 = "inflammatory response"

m1 = ["Tpi1","Pgk1","Pfkp","Pfkl","Hif1a","Gpi1","Gapdh","Aldoa",
      "Ldhb","Glut1","Eno1","Aco1","Aldoc","Slc2a1","Slc2a6","Hk1",
      "Hk2","Pfkp","Pgm1","Pgm2","Pkm2","Ldha"]  # 自身构造基因集合
m2 = ["Pparg", "Igfbp4", "Nos", "Cxcl1", "Jak2", "Trl6", "Myd88", "IL18", "Trl1", "Vcam1", "Icam1", "Ccl8",
      "Aim2", "Ccl2", "Ccl7", "Tnf", "Il6", "Il12a", "Il1b", "Sirt3"]  # 整体促进炎症
key_words = [key1, key2, key3, key4, key5, key6, key7, key8, key9, key10, key11, key12]
custom = m2
home_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
output_dir = os.path.join(home_dir, 'VOLCANO_HEATMAP_PREPROCESS')
input_name = 'GSE80185'
get_data(url, home_dir, input_name)
input_diff = os.path.join(home_dir, input_name+'.diff')
output_diff = os.path.join(output_dir, input_name+'.process.diff')

go_config = os.path.join(home_dir, 'config/gene2go.pickle')
with open(go_config, 'rb') as f:
    gene2go_dic = pickle.load(f)


# annotation
with open(output_diff, 'w') as fo:
    with open(input_diff) as f:
        i = 0
        scaned_set = set()
        # header_dic = {}
        # for idx, item in enumerate(header_line.split('\t')):
        #     header_dic[item] = idx
        for idx, line in enumerate(f):
            supp_list = [] # 追加各个基因是否属于某个信号通路
            flag = False
            line = line.replace('"', '')
            line_list = line.split('\t')[1:-1]

            if idx == 0:
                supp_list = ['custom'] + ['key'+str(i) for i in range(1, len(key_words)+1)]

            else:
                gene_symbol = line_list[-1]
                if not gene_symbol.strip():
                    continue

                if gene_symbol in custom:
                    supp_list.append('Y')
                else:
                    supp_list.append('N')

                if not gene_symbol in gene2go_dic:
                    go_infos = ('', '', '')
                else:
                    go_infos = gene2go_dic[gene_symbol]

                for key in key_words:
                    flag = 'N'
                    for go in go_infos:
                        if key in go:
                            flag = 'Y'
                            break
                    supp_list.append(flag)


                if gene_symbol in scaned_set:
                    continue

                scaned_set.add(gene_symbol)

            new_line = '\t'.join(line_list+supp_list)+'\n'
            fo.write(new_line)
    print(i)