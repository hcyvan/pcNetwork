with open('../reports/diff.lncRNA.646.list.csv') as f:
    rnas = [x.strip() for x in f.readlines()]
input_file='../../download-biomart/transcripts.unsplice.txt'
input_file='../../download-biomart/lncRNA646.transcript.unsplit.fasta'
with open(input_file) as f:
    p=False
    i=0
    for line in f:
        item=line.strip()
        if item.startswith('>'):
            #print(item)
            item_ids=item[1:].split('|')
            if item_ids[0] in rnas:
                i=i+1
                print(item)
                p=True
            else:
                p=False
        if p:
            continue
            print(item)




