# _*_ coding:utf-8 _*_

import os
import pickle
import numpy as np

# make anno dic
def build_id_convert_dic(inputs):
    with open(inputs) as f:
        has_headers = False
        id_convert_dic = {}
        for line in f:
            if line.startswith('#'):
                continue
            elif line.startswith('ID'):
                has_headers = True
                headers = line.strip().split('\t')
                header_dic = {}
                for idx, head in enumerate(headers):
                    header_dic[head] = idx
            elif has_headers:
                line_list = line.split('\t')
                gene_symbol = line_list[header_dic['Gene Symbol']]
                affy_id = line_list[header_dic['ID']]
                gobp = line_list[header_dic['Gene Ontology Biological Process']]
                gocc = line_list[header_dic['Gene Ontology Cellular Component']]
                gomf = line_list[header_dic['Gene Ontology Molecular Function']]
                id_convert_dic[affy_id] = [gene_symbol, gobp, gocc, gomf.strip()]
        return id_convert_dic




home_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
output_dir = os.path.join(home_dir, 'VOLCANO_HEATMAP_PREPROCESS')
input_name = 'GSE72518'
src = os.path.join(home_dir, input_name+'_series_matrix.txt')
tgt = os.path.join(output_dir, input_name+'.heatmap')

group_set = [0, 0, 0, 0, 0, 3, 3, 3, 1, 1, 1]  # GSE53986
group_name = ["ctrl", "ctrl", "ctrl", "ctrl", "lps", "lps", "lps", "lps"]


revise_group1 = [idx+1 for idx, item in enumerate(group_set) if item == 0]
revise_group2 = [idx+1 for idx, item in enumerate(group_set) if item == 1]
revise_group = revise_group1 + revise_group2


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
key12 = "oxidative phosphorylation"
key13 = "tricarboxylic acid cycle"

m1 = ["Tpi1","Pgk1","Pfkl","Hif1a","Gpi1","Gapdh","Aldoa",
      "Eno1","Slc2a1","Slc2a6","Hk1",
      "Hk2","Pgm2","Pkm2","Pkm","Ldha", "Mif",
      "Hk3", "Pgam1", "Fbp1", "Pgk2", "Pklr", "Ler3"]  # 自身构造基因集合

m2 = ["Nduf", "Core2", "Atp8", "Atp6", "Uqcrc1", "Cox4", "Cox5",
      "Cox7", "Atp5", "Atp7", "Ucp1", "Ak2", "Slc25a4", "Uqcc2", "Cox2"] #氧化磷酸化

key_list = [key1, key2, key3, key4, key5, key6, key7, key8, key9, key10, key11, key12, key13]  # 从GO描述中通过关键词筛选
custom_gene_list = m2

anno_cfg = os.path.join(home_dir, 'config/Affy_Mouse430_2.txt') # affy
anno_dic = build_id_convert_dic(anno_cfg)
# anno_cfg1 = os.path.join(home_dir, 'config/illum_mus_ref8_v2.pickle') # illu
# anno_cfg2 = os.path.join(home_dir, 'config/gene2go.pickle')
#
# with open(anno_cfg1, 'rb') as f:
#     anno_dic = pickle.load(f)
#
# with open(anno_cfg2, 'rb') as f:
#     go_dic = pickle.load(f)

with open(src) as f:
    with open(tgt, 'w') as fo:
        i = 0
        gene_set = set()
        header_line = ''
        sample_num = 0
        for line in f:
            if line.startswith('!') or not line.strip():
                continue
            line_list = line.strip().split('\t')
            probe_id = line_list[0].strip('"')
            line_list[0] = probe_id
            line_list[1:] = [line_list[idx] for idx in revise_group]
            if probe_id.startswith('ID'):
                line_list = line_list+["custom"]+["key"+str(i) for i in range(1, len(key_list)+1)]
                line_list[0] = "Gene"
                header_line = '\t'.join(line_list)+'\n'
                fo.write(header_line)
            else:
                supp_list = []
                anno_info = anno_dic.get(probe_id)
                if not anno_info:
                    print('failed')
                    i += 1
                    continue
                gene_symbol, gobp, gocc, gomf = anno_info
                if not gene_symbol.strip():
                    i += 1
                    continue

                gene_symbol = gene_symbol.split('/')[0].strip()
                if gene_symbol in gene_set:
                    print("failed2")
                    i += 1
                    continue

                gene_set.add(gene_symbol)

                data_list = np.array([float(x) for x in line_list[1:]])
                average = np.mean(data_list)
                norm_list = [str(x) for x in list(data_list/average)]

                for custom in custom_gene_list:
                    if gene_symbol.startswith(custom):
                        supp_list.append('Y')
                        break
                else:
                    supp_list.append('N')

                for key in key_list:
                    flag = 'N'
                    for go in [gobp, gocc, gomf]:
                        if key in go:
                            flag = 'Y'
                            break
                    supp_list.append(flag)
                line_list[0] = gene_symbol
                line_list = [line_list[0]] + norm_list + supp_list

                new_line = '\t'.join(line_list)+'\n'
                fo.write(new_line)
        print(i)



