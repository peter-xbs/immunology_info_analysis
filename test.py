# _*_ coding:utf-8 _*_

import numpy as np

genes = ["1460654_at","1460655_a_at",
         "1460656_a_at","1460657_at",
         "1460658_at","1460659_at",
         "1460660_x_at"]
# affy_dic_file = 'Affy_Mouse430_2.txt'
# gse_file = 'GSE53986_series_matrix.txt'
# gse_output = 'GSE53986_series_matrix_format.txt'

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


def convert_id_2_symbol_go(inputs, outputs, id_dic):
    with open(outputs, 'w') as fo:
        with open(inputs) as f:
            gene2go_dic = {}
            Set = set()
            for line in f:
                if line.startswith('!') or not line.strip():
                    continue
                line_list = line.strip().split('\t')
                affy_id = line_list[0].strip('"')
                if not affy_id:
                    continue

                if affy_id.startswith('ID'):
                    line_list[0] = 'Gene'
                    # line_list.insert(0, 'Gene')
                    # line_list.extend(['gobp', 'gocc', 'gomf'])
                    new_line = '\t'.join(line_list) + '\n'
                    fo.write(new_line)
                else:
                    if affy_id in id_dic:
                        gene_symbol, gobp, gocc, gomf = id_dic[affy_id]
                        gene2go_dic[gene_symbol] = (gobp, gocc, gomf)
                        gene_symbol = gene_symbol.split('///')[0].strip()
                        if gene_symbol in Set or not gene_symbol:
                            continue
                        Set.add(gene_symbol)
                        line_list[0] = gene_symbol
                        # line_list.insert(0, gene_symbol)
                        # line_list.extend([gobp, gocc, gomf])
                        new_line = '\t'.join(line_list) + '\n'
                        fo.write(new_line)
                    else:
                        print(affy_id)
            with open('gene2go.pickle', 'wb') as fo:
                import pickle
                pickle.dump(gene2go_dic, fo)

gse_output = 'GSE53986_series_matrix_format.txt'
gse_filt_output = 'GSE53986_series_matrix_format_filt.txt'
with open(gse_output) as f, open(gse_filt_output, 'w') as fo:
    # key_words = ['immun', 'inflam', 'b cell', 't cell', 'macrophage',
    #              'neutrophil', 'natural killer cell', 'dendritic cell',
    #              'mitochondrial', 'reactive oxygen species']
    key_words = ['macrophage', 'reactive oxygen species', 'mitochondrial']
    headers = [x.strip('"') for x in f.readline().split('\t')]
    header = headers[0]+'\t'+'\t'.join(headers[2:-3])+'\n'
    fo.write(header)
    gene_dic = {}
    for line in f:
        line_list = line.split('\t')
        gobp, gocc, gomf = line_list[-3], line_list[-2], line_list[-1]

        for key in key_words:
            if key in gobp or key in gocc or key in gomf:
                gene_symbol = line_list[0].split('///')[0]
                expr = np.array([float(x) for x in line_list[2:-3]])
                if gene_symbol not in gene_dic:
                    gene_dic[gene_symbol] = expr
                else:
                    gene_dic[gene_symbol] = (gene_dic[gene_symbol] + expr)/2
                new_line = line_list[0]+'\t'+'\t'.join(line_list[2:-3])+'\n'
        else:
            continue
    for gene in gene_dic:
        expr_s = '\t'.join([str(x) for x in gene_dic[gene]])
        new_line = '\t'.join([gene, expr_s])+'\n'
        fo.write(new_line)



if __name__ == '__main__':
    affy_dic_file = 'Affy_Mouse430_2.txt'
    gse_file = 'GSE53986_series_matrix.txt'
    gse_output = 'GSE53986_series_matrix_format2.txt'
    id_dic = build_id_convert_dic(affy_dic_file)
    convert_id_2_symbol_go(gse_file, gse_output, id_dic)





