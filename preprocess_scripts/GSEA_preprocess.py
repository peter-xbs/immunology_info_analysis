# _*_ coding:utf-8 _*_

import os

# 适用于2分组
output = '../GSEA_PREPROCESS'
if not os.path.exists(output):
    os.makedirs(output)

home_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
output_dir = os.path.join(home_dir, 'GSEA_PREPROCESS')
input_name = 'GSE72518'
src = os.path.join(home_dir, input_name+'_series_matrix.txt')
tgt = os.path.join(output_dir, input_name+'_series_matrix.gct')
group = os.path.join(output_dir,input_name+'_series_matrix.cls')

# group_set = [0,0,0,0,'x','x','x','x',1,1,1,1,'x','x','x','x'] # gse53986
# group_set = [0,0,0,1,1,1,'x','x','x','x','x','x'] # gse30552
# group_set = [0,0,0,0,0,1,1,1,1,1] # 27309z
# group_set = [1,1,1,0,0,0] # GSE43075
# group_name = ["LPS", "LPS", "LPS", "Control", "Control","Control"]
# group_comm = ["lps", "control"]
# group_set = [1, 1, 3, 3, 3, 3, 3, 0, 0, 0] # GSE60290
# group_name = ["lps", "lps", "control", "control", "control"]
# group_comm =["lps", "control"]
group_set = [0, 0, 0, 0, 0, 1, 1, 1, 3, 3, 3]  # GSE53986
group_name = ["ctrl", "ctrl", "ctrl", "ctrl", "ctrl", "IL4", "IL4", "IL4"]
group_comm = ["control", "lps"]

revise_group = [idx+1 for idx, item in enumerate(group_set) if item == 0 or item == 1]

with open(group, 'w') as fo:
    #group_control = [group_name[idx] for idx,item in enumerate(group_set) if item == 0]
    #group_exp = [group_name[idx] for idx, item in enumerate(group_set) if item == 1]
    #group_ = group_control + group_exp
    header_line = '\t'.join([str(len(group_name)), '2', '1'])+'\n'
    comm_line = '\t'.join(['#']+group_comm)+'\n'
    line = '\t'.join(group_name)+'\n'
    fo.write(header_line)
    fo.write(comm_line)
    fo.write(line)

with open(src) as f:
    with open(tgt, 'w') as fo:
        new_list = []
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
                sample_num = len(line_list)-1
                line_list = [item.strip('"') for item in line_list]
                line_list[0] = 'NAME'
                line_list.insert(1, 'Description')
                header_line = '\t'.join(line_list)+'\n'
            else:
                line_list.insert(1, 'NULL')
                new_line = '\t'.join(line_list)
                new_list.append(new_line)

        probe_set_num = len(new_list)
        first_line = '#1.2'+'\n'
        if sample_num == 0:
            print("Provide wrong input format!")
        second_line = '\t'.join([str(probe_set_num), str(sample_num)])+'\n'
        fo.write(first_line)
        fo.write(second_line)
        fo.write(header_line)
        for line in new_list:
            fo.write(line+'\n')


