import os


def handle_encode_tf_chip_seq_2015():
    data_file='data/ENCODE_TF_ChIP-seq_2015.txt'
    data_split=[['cell','tf','symbol'],]
    with open(data_file) as f:
        for item in f.readlines():
            data=item.strip().split('\t')
            if 'hg19' in data[0]:
                anno=data[0].split('_')
                data_split.extend([(anno[0],anno[1],x) for x in data[1:] if x != ''])
    for item in data_split:
        print('{}\t{}\t{}'.format(item[0], item[1], item[2]))


def handle_chea_2016():
    data_file='data/ChEA_2016.txt'
    data_split=[['tf','symbol'],]
    with open(data_file) as f:
        for item in f.readlines():
            data=item.strip().split('\t')
            if 'Human' in data[0]:
                anno=data[0].split('_')
                data_split.extend([(anno[0],x) for x in data[1:] if x != ''])
    for item in data_split:
        print('{}\t{}'.format(item[0], item[1]))

if __name__=='__main__':
    handle_chea_2016()
