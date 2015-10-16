#!/usr/bin/env python

import Image; 
import ImageDraw;
import ImageFont;
import string;
import math;
import sys;
import re;
import os;
usage="""
usage: drawcmp1allstate.py [-l LISTFILE] [-outpath DIR] [-q]
                           FILE [ FILE ... ]
Description: Draw general fractions of topology comparison classes
Options:
  -outpath DIR      set output path, (default: ./)
  -l  LISTFILE      set the fileListFile
  -pfamidlist FILE  Draw columns only for families within the supplied list file
  -colw  INT      set width for each column, (default: 1)
  -autosize y|n   whether use auto scale, (default: yes)
  -min-numseq INT set the minimal number of sequences for family 
                  to be selected, (default: 20)
  -max-numseq INT set the minimal number of sequences for family 
                  to be selected, (default: 10000)
  -min-numtm INT  set the minimal number of TMs for family 
                  to be selected, (default: 1)
  -max-numtm INT  set the minimal number of TMs for family 
                  to be selected, (default: 100)
  -addinfo        Add extra info to file name, $rootname.nseq$MIN-$MAX.nTM$MIN-$MAX.png

  -q              quiet mode
  -h, --help      print this help message and exit

Created 2011-09-16, updated 2011-11-25, Nanjiang Shu 
"""
#import comptopo as ct;
method_comparison=1;
colorsDict={}
#               IDT     OK     SHIFT  INV    INV_SHIFT    DIFF
colorsDict[1]=["green","red","pink","blue","lightgreen","black"];
#               IDT     SAME_NTERM   DIFF_NTERM
colorsDict[11]=["green","pink","blue"];

MIN_NUMSEQ = 20;
MAX_NUMSEQ = 10000;
MIN_NUMTM = 1;
MAX_NUMTM = 100;
isAddExtraInfoToFileName = False;

def PrintHelp():
    print usage;

def GetData(infile):#{{{
    data=[];      #  [ID, numseq, N1%, N2%, N3%, N4%, N5%, N6%]
    apd1=data.append;
    global method_comparison;
    try:
        fpin=open(infile,"r");
        lines=fpin.readlines();
        fpin.close()
        cnt=0;
        for line in lines:
            if len(line)>0 and line[0]!= '#':
                strs=line.split();
                if len(strs) == 17: #method_comparison == 1
                    method_comparison=1;
                    pfamid = strs[1]
                    numseq=int(strs[2]);
                    numTMIDT=int(strs[3]);
                    if (numseq >= MIN_NUMSEQ and numseq <= MAX_NUMSEQ and
                            numTMIDT >= MIN_NUMTM and numTMIDT <= MAX_NUMTM):
                        apd1([]);
                        apd2=data[cnt].append;
                        apd2(strs[1]);
                        apd2(int(strs[2]));
                        for i in range(11,17):
                            apd2(float(strs[i]));
                        cnt+=1;
                elif len(strs) == 11:
                    method_comparison=11;
                    numseq=int(strs[2]);
                    numTMIDT=int(strs[3]);
                    if (numseq >= MIN_NUMSEQ and numseq <= MAX_NUMSEQ and
                            numTMIDT >= MIN_NUMTM and numTMIDT <= MAX_NUMTM):
                        apd1([]);
                        apd2=data[cnt].append;
                        apd2(strs[1]);
                        apd2(int(strs[2]));
                        for i in range(8,11):
                            apd2(float(strs[i]));
                        cnt+=1;
                else:
                    sys.stderr.write("Unrecognized input datafile %s. Ignore.\n"
                            % infile);
                    break;


        return data
    except IOError:
        print >> sys.stderr, ("%s: Failed to read datafile %s. Ignore."
                % (sys.argv[0], infile));
    return data;
#}}}

def plotNumSeqDistribution(numSeqList, width, heightNumSeqPlotRegion,  #{{{
        marginleft, marginright, y, colwidth, draw):
    y0=y;
    maxNumSeq=max(numSeqList);
    numRecord=len(numSeqList);
    paddingtop=int(heightNumSeqPlotRegion*0.1+0.5);
    paddingbottom=int(heightNumSeqPlotRegion*0.05+0.5);
    x1=marginleft;
    y1=y0+paddingtop;
    x2=width-marginright;
    y2=y0+heightNumSeqPlotRegion-paddingbottom;
    box=[x1,y1,x2,y2];
    draw.rectangle(box, outline='black');

    # draw tics and text
    font_size=12;
    fnt = ImageFont.truetype(font_dir+font, font_size)
    step=max(1,maxNumSeq/3);
    ss=str(step);
    step=int(ss[0]);
    for i in range(len(ss)-1):
        step*=10;
    ytic=0;
    lengthtic=int(width*0.01+0.5);
#    step=200;
    while ytic < maxNumSeq:
        x1=marginleft-lengthtic;
        y1=(y0 + heightNumSeqPlotRegion-paddingbottom - int(ytic/float(maxNumSeq)
            * (heightNumSeqPlotRegion - paddingbottom-paddingtop)+0.5));
        x2 = x1 + lengthtic;
        y2=y1;
        draw.line([x1, y1, x2, y2],fill="black");
        text = "%d"%ytic;
        (textWidth,textHeight) = fnt.getsize(text);
        draw.text((x1-textWidth-lengthtic-3,y1-textHeight/2), text, font=fnt,
                fill='black');
        ytic+=step;
    

    hList = [int(n/float(maxNumSeq) *
        (heightNumSeqPlotRegion-paddingbottom-paddingtop) + 0.5) for n in
        numSeqList];
    x=marginleft;
    sizeSquare=4
    pointList=[];
    papd=pointList.append;
    for h in hList:
        x1=x+colwidth/2-sizeSquare/2;
        y1=y0+heightNumSeqPlotRegion-paddingbottom-h-sizeSquare/2;
        x2=x1+sizeSquare;
        y2=y1+sizeSquare;
        box=[x1,y1,x2,y2];
        draw.ellipse(box, outline='red');
        papd((x+colwidth/2,y0+heightNumSeqPlotRegion-paddingbottom-h));
        x+=colwidth;
    for i in xrange(0,len(pointList)-1,1):
        draw.line([pointList[i],pointList[i+1]],fill="green");

#}}}

def AutoSize(numRecord):#{{{
    marginleft=5;
    marginright=5;
    margintop=5;
    marginbottom=30;
    colwidth=1;
    heightDrawRegion=1000;
    ratioHeight2Width=0.618;
    minSize=500;
    minHeight=int(minSize*ratioHeight2Width+0.5);
    font_size=14;

    if numRecord > minSize:
        colwidth=1;
        heightDrawRegion = int(colwidth*numRecord*ratioHeight2Width+0.5);
    else:
        colwidth=int(math.ceil(minSize/float(numRecord)));
        heightDrawRegion = int (colwidth*numRecord*ratioHeight2Width+0.5);
    marginleft   = max(1,int(colwidth*numRecord*0.12+0.5));
    marginright  = max(1,int(colwidth*numRecord*0.05+0.5));
    margintop    = max(1,int(heightDrawRegion*0.01+0.5));
    marginbottom = max(25,int(heightDrawRegion*0.05+0.5));
    font_size+=int((heightDrawRegion/minHeight-1)*1.5) ;
    fnt = ImageFont.truetype(font_dir+font, font_size);

    return (colwidth,heightDrawRegion, marginleft, marginright,margintop,marginbottom, font_size, fnt);
#}}}
def plotcolumn(values, x,y, heightDrawRegion, colwidth, draw):#{{{
    pixheights=[int(j*heightDrawRegion/100+0.5) for j in values];
    #print pixheights;
    numBin=len(pixheights);
    y0=y;   #preserve the initial y
    x1=x+colwidth;
    for i in range(numBin):
        h=pixheights[i];
        y1=y+h;
        if i==numBin-1:
            y1=heightDrawRegion+y0;
        box=[x,y,x1,y1];
        draw.rectangle(box, fill=colorsDict[method_comparison][i]);
        y+=h;
#}}}
def DrawCMPAllStateImage(infile, font_size, fnt):#{{{
    marginleft=5;
    marginright=5;
    margintop=5;
    marginbottom=30;
    colwidth=1;
    heightDrawRegion=1000;

    rootname=os.path.basename(os.path.splitext(infile)[0]);
    outfile="";
    if isAddExtraInfoToFileName:
        outfile= (outpath+os.sep+"%s.nseq%d-%d.nTM%d-%d.png"%(rootname,
            MIN_NUMSEQ, MAX_NUMSEQ, MIN_NUMTM, MAX_NUMTM));
    else:
        outfile=outpath+os.sep+"%s.png"%(rootname);

    data = GetData(infile);
    numRecord=len(data);
    if numRecord > 0:
        if isAutoSize:
            (colwidth, heightDrawRegion, marginleft, marginright,
                    margintop,marginbottom, font_size, fnt) = AutoSize(numRecord);

        heightNumSeqPlotRegion=heightDrawRegion/3;
        width=numRecord*colwidth+marginleft+marginright;
        height = (heightDrawRegion + heightNumSeqPlotRegion + margintop +
                marginbottom);
#        mode="P";
        mode="RGB";

#        img=Image.new(mode, (width, height),255);

# draw rotated text
        rotated_img=Image.new(mode, (height, width),"white");
        #rotated_img=img.rotate(-90, expand=True);
        draw = ImageDraw.Draw(rotated_img); # setup to draw on the rotated img
        text="Number of sequences";
        fnt = ImageFont.truetype(font_dir+font, 12);
        (textWidth,textHeight) = fnt.getsize(text);
        x=marginbottom + heightNumSeqPlotRegion/2 - textWidth/2;
        y = 5;
        draw.text((x,y), text, font=fnt, fill='black');
        img=rotated_img.rotate(90, expand=True);

#        img.show();

        draw1 = ImageDraw.Draw(img); # setup to draw on the main image

        x=marginleft;
        y=margintop;
        for item in data:
            #print item[1:];
            plotcolumn(item[2:], x,y, heightDrawRegion, colwidth, draw1);
            x+=colwidth;
            y=margintop;

# plot numseq distributions
        numSeqList=[];
        apd=numSeqList.append;
        for itm in data:
            apd(itm[1]);
        x=marginleft;
        y=margintop + heightDrawRegion;
        plotNumSeqDistribution(numSeqList, width, heightNumSeqPlotRegion,
                marginleft, marginright, y, colwidth, draw1)

        text="%d families"%numRecord;
        (textWidth,textHeight) = fnt.getsize(text);
        x = (width-textWidth)/2;
        y = (heightDrawRegion + heightNumSeqPlotRegion +
                margintop+marginbottom/2 - textHeight/2);
        fg='black';
        draw1.text((x,y), text, font=fnt, fill=fg);

        #img.show();
        img.save(outfile);
        if not isQuiet:
            print outfile, "has been output";
#}}}
if __name__ == "__main__":

    numArgv=len(sys.argv)
    if numArgv < 2:
        PrintHelp()
        sys.exit()

    isQuiet=False;
    isAutoSize=True;
    outpath="./";
    fileListFile=None;
    fileList=[];
    #infile="result_mode2_uniq/all.sorted.generalinfo.new.txt";


    i = 1;
    isNonOptionArg=False
    while i < numArgv:#{{{
        if isNonOptionArg == True:
            fileList.append(sys.argv[i]);
            isNonOptionArg=False;
            i += 1;
        elif sys.argv[i] == "--":
            isNonOptionArg=True;
            i += 1;
        elif sys.argv[i][0] == "-":
            if sys.argv[i] ==  "-h" or  sys.argv[i] == "--help":
                PrintHelp();
                sys.exit();
            elif sys.argv[i] == "-outpath" or sys.argv[i] == "--outpath":
                outpath=sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-colw" or sys.argv[i] == "--colw":
                colwidth=int(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] == "-l" :
                fileListFile=sys.argv[i+1];
                i += 2;
            elif sys.argv[i] == "-min-numseq" :
                MIN_NUMSEQ=int(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] == "-max-numseq" :
                MAX_NUMSEQ=int(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] == "-min-numtm" :
                MIN_NUMTM=int(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] == "-max-numtm" :
                MAX_NUMTM=int(sys.argv[i+1]);
                i += 2;
            elif sys.argv[i] == "-autosize" or sys.argv[i] == "--autosize" :
                if sys.argv[i+1].lower()[0]== "y":
                    isAutoSize=True;
                else:
                    isAutoSize=False;
                i += 2;
            elif sys.argv[i] == "-q":
                isQuiet=True;
                i += 1;
            elif sys.argv[i] == "-addinfo":
                isAddExtraInfoToFileName=True;
                i += 1;
            else:
                print >> sys.stderr, "Error! Wrong argument:", sys.argv[i];
                sys.exit()
        else:
            fileList.append(sys.argv[i])
            i += 1
#}}}
    if fileListFile != None:
        try:
            fp=open(fileListFile,"r");
            fileList+=fp.read().split();
            fp.close();
        except IOError:        
            print >> sys.stderr, "file %s does not exist."

    font_dir= "/usr/share/fonts/truetype/ttf-dejavu/";
    font="DejaVuSansMono.ttf";
    font_size = 12;
    fnt = ImageFont.truetype(font_dir+font, font_size)
    os.system("mkdir -p %s"%(outpath));

    for infile in fileList:
        DrawCMPAllStateImage(infile, font_size, fnt);
    
