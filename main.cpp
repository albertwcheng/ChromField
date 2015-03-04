/***************************************************************************
 Copyright 2010 Wu Albert Cheng <albertwcheng@gmail.com>
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 *******************************************************************************/ 



#include <iostream>
#include <fstream>
#include <limits.h>
#include <vector>
#include <deque>
#include <map>
#include <set>
#include <list>
#include <queue>
#include <algorithm>
#include <BitString.h>
#include <StringUtil.h>
#include "ChromField.h"





void printFGHelp(const char* programname)
{
	cerr<<"Usage:"<<programname<<" command (see below)"<<endl;
	cerr<<"Description: Operates on ChromField files"<<endl;
	cerr<<"Programs:"<<endl;
	
	cerr<<"BedToChromField <chrSizeList> <outChromField> <inbed1> <inbed2> ... <inbedN>"<<endl;
	cerr<<"\tGenerete binary ChromField file [0 for non-interval, 1 for interval]"<<endl;
    
	cerr<<"cBedToChromField <chrSizeList> <outChromField> <inbed1> <inbed2> ... <inbedN>"<<endl;
	cerr<<"\tGenerete binary complement ChromField file  [1 for non-interval, 0 for interval]"<<endl;
    
    cerr<<"ChromFieldToBed <chrSizeList> <inChromField> <outBed>"<<endl;
    cerr<<"\tExtract bed intervals from ChromField file"<<endl;
    
	cerr<<"SelectBedItemsByOverlap <chrSizeList> <inbed> <inChromField> <overlapSize; -1 means bed item is contained; 0 means get everything> <outBed> [attachOverlapLen attachBitString[S]  ...]"<<endl;
	cerr<<"\tOutput bed files for reads occuring [thresholdLow,thresholdHigh] inclusive times "<<endl;
    
    cerr<<"Union <chrSizeList> <outChromField> <inChromField1> <inChromField2> ... <inChromFieldN>"<<endl;
    cerr<<"\tOutput the union of inChromField files"<<endl;
    
    cerr<<"Intersect <chrSizeList> <outChromField> <inChromField1> <inChromField2> ... <inChromFieldN>"<<endl;
    cerr<<"\tOutput the intersection of inChromField files"<<endl;
    
    cerr<<"Subtract <chrSizeList> <outChromField> <inChromField1> <inChromField2> ... <inChromFieldN>"<<endl;
    cerr<<"\tSubtract inChromField1 by inChromField2...N and outut to outChromField file"<<endl;
    

	//cerr<<"-pk binary"<<endl;
	//cerr<<"\tPrint the content of a binary KEYEDPOSITION file"<<endl;
	
}


string reverseString(const string& src){
    string reta;
    for(int i=src.length()-1;i>=0;i--){
        reta+=src[i];
    }
    
    return reta;
}

void SelectBedItemsByOverlap(int argc,const char**argv){
    //SelectBedItemsByOverlap <chrSizeList> <inbed> <inChromField> <overlapSize; -1 means bed item is contained; 0 means get everything> <outBed> [ attachOverlapLen attachBitString[S] ...]
	
    if(argc<7)
	{
		printFGHelp(argv[0]);
		return;
	}
    
   
    
    int64_t overlapSizeLower=StringUtil::atoi(argv[5]);
    
    ifstream fil(argv[3]);
    string line;
    vector<string> splits;
    
    bool attachBitString=false;
    bool attachOverlapLen=false;
    bool strandSpec=false;
    
    for(int c=7;c<argc;c++){
        string sargv(argv[c]);
        if(sargv=="attachOverlapLen"){
            attachOverlapLen=true;
        }else if(sargv=="attachBitString"){
            attachBitString=true;
        }else if(sargv=="attachBitStringS"){
            attachBitString=true;
            strandSpec=true;
        }
    }
    
    ChromField chromfield(argv[4],argv[2]);
    
    ofstream ofil(argv[6]);
    uint64_t lino=0;
    
    while(fil.good())
    {
        lino+=1;
        if(lino%1000000==1){
            cerr<<"processing inbed line "<<lino<<endl;
        }
        getline(fil,line);
        if(fil.good() && line.length()>0){
            
            //cerr<<"a"<<endl;
            StringUtil::split(line,"\t",splits);
            if(splits.size()<3){
                continue;
            }
            string chrom=splits[0];
            uint64_t start0=StringUtil::atoi(splits[1]);
            uint64_t end1=StringUtil::atoi(splits[2]);
            
            /*map<string,chromcoord>::iterator chromcoordI=chromMap.find(chrom);
            if(chromcoordI==chromMap.end()){
                cerr<<"Error: Chrom "<<chrom<<" not defined in chromSizes file. line skipped: "<<line<<endl;
                continue;
            }
            if(start0<0){
                cerr<<"Error: bed entry start0<0 thus invalid. line skipped: "<<line<<endl;
                continue;
            }
            if(start0>=end1){
                cerr<<"Error: bed entry start0>=end1 thus invalid. line skipped: "<<line<<endl;
                continue;
            }
            
            if(end1>chromcoordI->second.length()){
                cerr<<"Error: bed entry has end1 > length of chrom. line skipped: "<<line<<endl;
                continue;
            }
            
            uint64_t chromBitStart0=chromcoordI->second.start0;*/
           //  cerr<<"b"<<endl;
            //uint64_t olen;
            
            
            
            pair<string,uint64_t> bsOL=chromfield.getBits(chrom,start0,end1);
            uint64_t olen=bsOL.second;
            
            //=chromfield.getOverlapLength(chrom,start0,end1);
          
            
            //   cerr<<"c"<<endl;
            if(overlapSizeLower==-1 && olen<end1-start0){
                continue; //need complete containment;
            }
            
            if(olen<overlapSizeLower){
                continue; //not met;
            }
            
            
            ofil<<line;
            
            if(attachOverlapLen){
                ofil<<"\t"<<olen;
            }
            
            if(attachBitString){
                if(strandSpec){
                    string& strand=splits[5];
                    if(strand=="+"){
                        ofil<<"\t"<<bsOL.first;
                    }else{
                        ofil<<"\t"<<reverseString(bsOL.first);
                    }
                }else{
                    ofil<<"\t"<<bsOL.first;
                }
            }
            
            ofil<<endl;
        }
    }
    
    ofil.close();

}
void BedToChromField(int argc,const char** argv)
{
	if(argc<5)
	{
		printFGHelp(argv[0]);
		return;
	}

	//setBitOverBedIntervals <chrSizeList> <outChromField>  <inbed1> <inbed2> ... <inbedN>
    ChromField chromfield(argv[2],0);
    
    //bool value=(string(argv[4])!="0");
    
    for(int i=4;i<argc;i++){
        cerr<<"processing bedfile "<<argv[i]<<endl;
        chromfield.setBitsOnBedIntervals(argv[i],true);
    }
    
    //chromfield.printBed(cout);
    
    chromfield.writeToFile(argv[3]);
}

void cBedToChromField(int argc,const char** argv)
{
	if(argc<5)
	{
		printFGHelp(argv[0]);
		return;
	}
    
	//setBitOverBedIntervals <chrSizeList> <outChromField>  <inbed1> <inbed2> ... <inbedN>
    ChromField chromfield(argv[2],255); //11111111
    
    //bool value=(string(argv[4])!="0");
    
    for(int i=4;i<argc;i++){
        cerr<<"processing bedfile "<<argv[i]<<endl;
        chromfield.setBitsOnBedIntervals(argv[i],false);
    }
    
    //chromfield.printBed(cout);
    
    chromfield.writeToFile(argv[3]);
}

void IntersectChromFields(int argc,const char** argv)
{
    if(argc<6)
    {
        printFGHelp(argv[0]);
        return;
    }
    
    //Intersect <chrSizeList> <outChromField> <inChromField1> <inChromField2> ... <inChromFieldN>
    ChromField leftChromField(argv[4],argv[2]); //load the first chrom field from file
    
    
    for(int i=5;i<argc;i++){  //starting from inChromFile2
        cerr<<"processing chromField file "<<argv[i]<<endl;
        ChromField *rightChromField=new ChromField(argv[i],argv[2]);
        leftChromField&=(*rightChromField);
        delete rightChromField;
    }
    
    leftChromField.writeToFile(argv[3]);
}


void UnionChromFields(int argc,const char** argv)
{
    if(argc<6)
    {
        printFGHelp(argv[0]);
        return;
    }
    
    //Union <chrSizeList> <outChromField> <inChromField1> <inChromField2> ... <inChromFieldN>
    ChromField leftChromField(argv[4],argv[2]); //load the first chrom field from file
    
    
    for(int i=5;i<argc;i++){  //starting from inChromFile2
        cerr<<"processing chromField file "<<argv[i]<<endl;
        ChromField *rightChromField=new ChromField(argv[i],argv[2]);
        leftChromField|=(*rightChromField);
        delete rightChromField;
    }
    
    leftChromField.writeToFile(argv[3]);
}

void SubtractChromFields(int argc,const char** argv)
{
    if(argc<6)
    {
        printFGHelp(argv[0]);
        return;
    }
    
    //Subtract <chrSizeList> <outChromField> <inChromField1> <inChromField2> ... <inChromFieldN>
    ChromField leftChromField(argv[4],argv[2]); //load the first chrom field from file
    
    
    for(int i=5;i<argc;i++){  //starting from inChromFile2
        cerr<<"processing chromField file "<<argv[i]<<endl;
        ChromField *rightChromField=new ChromField(argv[i],argv[2]);
        leftChromField-=(*rightChromField);
        delete rightChromField;
    }
    
    leftChromField.writeToFile(argv[3]);
}


void ChromFieldToBed(int argc,const char** argv)
{
	if(argc<5)
	{
		printFGHelp(argv[0]);
		return;
	}
    
	//chromFieldToBed <chrSizeList> <inChromField> <outBed>
    ChromField chromfield(argv[3],argv[2]);
    ofstream outBed(argv[4]);
    chromfield.printBed(outBed);
    outBed.close();
 }

int main(int argc, const char **argv)
{

	cerr<<"ChromField"<<endl;
	cerr<<"[Built:"<<__DATE__<<" "<<__TIME__<<"]"<<endl;
	if(argc<2 || !strcmp(argv[1],"-h"))
	{
		printFGHelp(argv[0]);
	}
    else if(!strcmp(argv[1],"BedToChromField"))
	{
        BedToChromField(argc,argv);
        
	}else if(!strcmp(argv[1],"cBedToChromField"))
	{
        cBedToChromField(argc,argv);
        
	}
    else if(!strcmp(argv[1],"ChromFieldToBed"))
    {
        ChromFieldToBed(argc,argv);
    }
    else if(!strcmp(argv[1],"SelectBedItemsByOverlap")){
        SelectBedItemsByOverlap(argc,argv);
    }else if(!strcmp(argv[1],"Intersect")){
        IntersectChromFields(argc,argv);
    }else if(!strcmp(argv[1],"Union")){
        UnionChromFields(argc,argv);
    }else if(!strcmp(argv[1],"Subtract")){
        SubtractChromFields(argc,argv);
    }

    
    cerr<<"<Done>"<<endl;
	return 0;

}


