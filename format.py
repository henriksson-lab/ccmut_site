color_space='gray';
color_orig='blue';
color_rep='red';
color_neither='black';

def getcolor(c, orig, rep):
    if c=='-':
        return color_space
    elif c==orig:
        return color_orig
    elif c==rep:
        return color_rep
    else:
        return color_neither

def color_char(col, c):
    if col=='':
        return c
    else:
        return '</span><span style=\"color:'+col+';\">'+c

def redcolor(before,after):
    if before==after:
        return ''
    else:
        return after

def color_oneseq(oneseq,origseq,repseq):
    #Colors to use
    colors_for_seq=[getcolor(x[0],x[1],x[2]) for x in zip(list(oneseq), origseq, repseq)]

    #Figure out which color transitions are not needed - saves bandwidth
    colors_before = colors_for_seq[0:(len(colors_for_seq)-1)]
    colors_before.insert(0,"black")
    colors_red = [redcolor(x[0],x[1]) for x in zip(colors_before, colors_for_seq)]

    #Color the sequence
    colored_seq = ''.join([color_char(x[0],x[1]) for x in zip(colors_red,list(oneseq))])
    return '<span style=\"color:black;\">' + colored_seq + '</span>'


#######################


#lines=lines[0:3]
#origseq=list(lines[2])
#repseq=list(lines[1])
#oneseq=lines[0]
#co = color_oneseq(oneseq,origseq,repseq)
#print(co)

def color_alignment(lines):
    lines=lines.rstrip().split("\n")

    repseq   = list(lines[len(lines)-1])
    origseq  = list(lines[len(lines)-2])

    colseqs = [color_oneseq(oneseq,origseq,repseq) for oneseq in lines]

    return '<br/>'.join(colseqs)



def main():
    f = open("ex.fasta", "r")
    lines=f.read()
    print(color_alignment(lines))

if __name__ == "__main__":
    main()
