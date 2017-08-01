import sys

if len(sys.argv)<4:
    print("usage: fastmd.py reference.fasta input.sam output.sam")
    exit()

def readref(filename):
    ref={}
    seqid=None
    for line in open(filename,"r").readlines():
        line=line.strip()
        if line=="":
            continue
        elif line[0]==">":
            seqid=line[1:].split(" ")[0].split("\t")[0]
            ref[seqid]=[]
        else:
            ref[seqid].append(line.upper())

    for seqid in ref:
        ref[seqid]="".join(ref[seqid])
    return ref

def cigar2ops(cigar):
    ops=[]
    buf=""
    for char in cigar:
        if char.isdigit():
            buf+=char
        else:
            ops.append((char,int(buf)))
            buf=""
    if len(buf):
        raise Exception("Leftover position in CIGAR")
    return ops

def reconstructalignment(pos,cigar,ref,read):
    refaln=""
    readaln=""
    refpos=pos
    readpos=0
    for op,length in cigar2ops(cigar):
        if op=="M" or op=="=" or op=="X":
            refaln+=ref[refpos:refpos+length]
            refpos+=length
            readaln+=read[readpos:readpos+length]
            readpos+=length
        elif op=="I":
            refaln+="-"*length
            readaln+=read[readpos:readpos+length]
            readpos+=length
        elif op=="D" or op=="N":
            refaln+=ref[refpos:refpos+length]
            refpos+=length
            readaln+="-"*length
        elif op=="S":
            readpos+=length
        elif op=="H" or op=="P":
            pass
        else:
            raise Exception("Unhandled CIGAR operation: %s"%op)
    if len(readaln)!=len(refaln):
        raise Exception("Alignment reconstruction failed")
    return refaln,readaln

def alignment2md(refaln,readaln):
    md=""
    counter=0
    starteddel=False
    for i in range(len(refaln)):
        if refaln[i]==readaln[i]:
            counter+=1
            starteddel=False
        elif refaln[i]=="-":
            starteddel=False
        else:
            if readaln[i]=="-":
                if not starteddel:
                    md+=str(counter)
                    md+="^"
                    starteddel=True
                md+=refaln[i]
                counter=0
            else:
                md+=str(counter)
                md+=refaln[i]
                counter=0
                starteddel=False
    return md+str(counter)

def verify(tags,md):
    origmd=None
    for tag in tags:
        if "MD:Z:" in tag:
            origmd=tag[5:]
    return md==origmd

if __name__ == "__main__":
    ref=readref(sys.argv[1])
    outh=open(sys.argv[3],"w")
    lid=0

    for line in open(sys.argv[2],"r"):
        lid+=1
        if line.strip()=="" or line[0]=="@":
            outh.write(line)
            continue
        parts=line.split("\t")
        if len(parts)<11:
            print("Warning: Line %d in input is invalid"%lid)
            continue
        qname,flag,rname,pos,mapq,cigar,rnext,pnext,tlen,readseq,qual=parts[:11]
        tags=parts[11:]
        readseq=readseq.upper()
        if not rname in ref:
            print("Warning: Missing sequence '%s' in reference for read %s"%(rname,qname) )
            continue
        refseq=ref[rname]

        try:
            refaln,readaln=reconstructalignment(int(pos)-1,cigar,refseq,readseq)
        except Exception as e:
            print("Warning: Alignment reconstruction failed for read %s: %s"%(qname,str(e)))
            continue

        try:
            md=alignment2md(refaln,readaln)
        except Exception as e:
            print("Warning: MD string computation failed for read %s: %s"%(qname,str(e)))
            continue

        outh.write(line.strip()+"\tMD:Z:"+md+"\n")

    outh.close()
