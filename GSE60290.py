# _*_ coding:utf-8 _*_

import requests

#url = 'https://www.ncbi.nlm.nih.gov/geo/geo2r/backend/geo2r.cgi?ctg_time=1538533969&job_key=JSID_01_588503_130.14.22.10_9000_geo2r'
# url = 'https://www.ncbi.nlm.nih.gov/geo/geo2r/backend/geo2r.cgi?ctg_time=1540824113&job_key=JSID_01_755259_130.14.18.128_9000_geo2r'
# data = requests.get(url).content.decode()
#
# with open('gse60290.txt', 'w') as fo:
#     fo.write(data)
gene_set = set()
with open('GSE60290', 'w') as fo:
    with open('gse60290.txt') as f:
        i = 0
        j = 0
        for line in f:
            line_list = line.strip('\n').split('\t')
            line_list = line_list[1:7]
            line_list = [item.strip('"') for item in line_list]
            gene = line_list[-1].strip()
            if not gene:
                j += 1
                new_line = '\t'.join(line_list)+'\n'
                fo.write(new_line)
            elif gene not in gene_set:
                gene_set.add(gene)
                new_line = '\t'.join(line_list)+'\n'
                fo.write(new_line)
            else:
                i += 1
        print(i,j)
#
# with open('GSE53986_series_matrix_format_filt.txt') as f:
#     gene_set = set()
#     for line in f:
#         line_list = line.split('\t')
#         gene = line_list[0]
#         gene_set.add(gene)
#
# with open('gse60290.txt') as f:
#     for line in f:
#         print(line)
#
# from time import sleep, ctime
# import threading
#
# loops = [4,2]
#
# def loop(nloop, nsec):
#     print('starting at:', nloop, 'at:', ctime())
#     sleep(nsec)
#     print('loop', nloop, 'done at: ', ctime())
#
# def main():
#     print('starting at:', ctime())
#     threads = []
#     nloops = range(len(loops))
#
#     for i in nloops:
#         t = threading.Thread(target=loop, args=(i, loops[i]))
#         threads.append(t)
#
#     for i in nloops:
#          threads[i].start()
#     for i in nloops:
#         threads[i].join()
#
#     print('all done at:', ctime())
#
#
# if __name__ == '__main__':
#     main()
#


