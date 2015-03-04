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
 
#ifndef __CHROM_FIELD_H
#define __CHROM_FIELD_H

#include <utility>
#include <BitString.h>
#include <stdint.h>
#include <StringUtil.h>
#include <map>
#include <vector>
#include <string>
using namespace std;

class ChromField:public BitString
{
	
    
public:
    
    class chromcoord{
    public:
        string chrom;
        uint64_t start0;
        uint64_t end1;
        chromcoord():start0(-1),end1(-1){}
        chromcoord(const string& _chrom,uint64_t _start0,uint64_t _end1):chrom(_chrom),start0(_start0),end1(_end1){}
        uint64_t length(){
            return end1-start0;
        }
    };
    
	map<string,chromcoord> chromMap;
    vector<chromcoord> orderedChromBitRange;
    
	uint64_t totalLength;
	
	uint64_t readChromMap(const string& chromSizeFile)
	{
		ifstream fil(chromSizeFile.c_str());
		string line;
		vector<string> splits;
		uint64_t start0=0;
		
        cerr<<"read chromSizeFile "<<chromSizeFile<<endl;
        
		while(fil.good()){
			line="";
			getline(fil,line);
			if(fil.good() && line.length()>0){
                
				StringUtil::split(line,"\t",splits);
				//cerr<<line<<"\t"<<splits.size()<<endl;
                if(splits.size()<2)
					continue;
				
				
				string chrom=string(splits[0]);
				uint64_t len=StringUtil::atoi(splits[1]);
                
                cerr<<chrom<<" len="<<len<<endl;
                chromcoord cc(chrom,start0,start0+len);
                
                
				chromMap.insert(map<string,chromcoord>::value_type(chrom,cc));
                //cerr<<chrom<<" blen="<<len<<endl;

                orderedChromBitRange.push_back(cc);
				start0+=len;	
				//cerr<<chrom<<" aaaa="<<len<<endl;
			}
		}
		
		totalLength=start0;
        //cerr<<"totalLength="<<totalLength<<endl;
		return totalLength;		
	}
	ChromField(const string& chromSizeFile,Byte valuePerByte)
	{
        this->init(readChromMap(chromSizeFile),valuePerByte);
		//this->BitString::print(cerr,true);
	}
	ChromField(const string& chromFieldFile, const string&chromSizeFile):BitString(chromFieldFile){
		readChromMap(chromSizeFile);
	}
	
	ChromField(uint64_t _totalLength):BitString(_totalLength),totalLength(_totalLength)
	{
		
	}
	
	/*bool getBit(const string& chrom,uint64_t coord0)
	{
        uint64_t chromBitStart0=chromMap[chrom].start0;
		return this->BitString::getBit(chromBitStart0+coord0);
	}*/
	
    void setBitsOnBedIntervals(const string& filename,bool value){
        ifstream fil(filename.c_str());
        string line;
        vector<string> splits;
        
       
        
        while(fil.good())
        {
            getline(fil,line);
            if(fil.good() && line.length()>0){
                StringUtil::split(line,"\t",splits);
                string chrom=splits[0];
                uint64_t start0=StringUtil::atoi(splits[1]);
                if(splits[1][0]=='-'){
                    start0=0;
                    cerr<<"Warning: Start0<0, trimed at 0. line: "<<line<<endl;
                }
                uint64_t end1=StringUtil::atoi(splits[2]);
                
                map<string,chromcoord>::iterator chromcoordI=chromMap.find(chrom);
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
                    //cerr<<"Error: bed entry has end1 > length of chrom. line skipped: "<<line<<endl;
                    cerr<<"Warning: bed entry has end1 > length of chrom. trimed to chrom length="<<chromcoordI->second.length()<<" line: "<<line<<endl;
                    end1=chromcoordI->second.length();
                    //continue;
                }
                
                uint64_t chromBitStart0=chromcoordI->second.start0;
                
                this->setBitsOnRange(chromBitStart0+start0,chromBitStart0+end1,value);
            }
        }
        
        //cerr<<BitString::numBytes<<endl;
        //this->BitString::print(cout,true);
        
        fil.close();
        
     }
    
    uint64_t getOverlapLength(const string& chrom,uint64_t start0,uint64_t end1){
        uint64_t olen=0;
        map<string,chromcoord>::iterator chromcoordI=chromMap.find(chrom);
        if(chromcoordI==chromMap.end()){
            cerr<<"Error: Chrom "<<chrom<<" not defined in chromSizes file. getOverlapLength operation aborted"<<endl;
            return 0;
        }
        
        uint64_t chrom0=chromcoordI->second.start0;
        for(uint64_t i=chrom0+start0;i<chrom0+end1;i++)
        {
            olen+=getBit(i);
        }
        return olen;
    }
    
    pair<string,uint64_t> getBits(const string& chrom,uint64_t start0,uint64_t end1){
        pair<string,uint64_t> bsOL;
        bsOL.second=0;
        uint64_t olen=0;
        
        map<string,chromcoord>::iterator chromcoordI=chromMap.find(chrom);
        if(chromcoordI==chromMap.end()){
            cerr<<"Error: Chrom "<<chrom<<" not defined in chromSizes file. getOverlapLength operation aborted"<<endl;
            //return 0;
            return bsOL;
        }
        
        uint64_t chrom0=chromcoordI->second.start0;
        for(uint64_t i=chrom0+start0;i<chrom0+end1;i++)
        {
            bool thisBit=getBit(i);
            if(thisBit){
                olen++;
                bsOL.first+="1";
            }else{
                bsOL.first+="0";
            }
           
            
        }
        
        
        bsOL.second=olen;
        return bsOL;
        
    }
    
    
    void printBed(ostream& os)
	{
        

        
        for(vector<chromcoord>::iterator B=orderedChromBitRange.begin();B!=orderedChromBitRange.end();B++){
          
            
            cerr<<"printing chrom "<<B->chrom<<endl;
            uint64_t start0=0;
            uint64_t end1=0;
            for(uint64_t i=B->start0;i<B->end1;i++){
                bool value=this->BitString::getBit(i);

                if(value){ //"1"
                    if(end1==0){
                        start0=i;
                    }
                    end1=i+1;
                }else  if(end1>0){ //"0" and previously "1"
                    //out
                    os<<B->chrom<<"\t"<<(start0-B->start0)<<"\t"<<(end1-B->start0)<<endl;
                    end1=0;
                }
            }
            
            if(end1>0){ //final print for this chromosome
                os<<B->chrom<<"\t"<<(start0-B->start0)<<"\t"<<(end1-B->start0)<<endl;
            }
            
        }
    }
	
};
 
 
 
#endif