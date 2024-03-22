///说明：seq用2bit存，Occ用popCnt找
///修改：2023.8.28，基于glocate12，将bwtcpy放入occ中，调优
///修改：2023.9.6，去掉STL中的map，改用<开始排位，个数>，线性扫描
///修改：2023.11.5，区间排序，提高mp的效率
///修改：2024.3.9，保存找到的read及seeds信息，方便复现
///修改：2024.3.12，读入已找到的read集合文件
///修改：2024.3.17，将C表，pOcc和sSA的long改为uint
///修改：2024.3.17，将pOcc和cBWT交叉存
///修改：2024.3.19，进一步对Effaln优化
///修改：2024.3.22，准备git发布版本
#include <bits/stdc++.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string>
#include <math.h>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <bitset>
#include <vector>
#include <unistd.h>
#include <algorithm>
#include <unistd.h>
#include <sys/mman.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include "divsufsort/divsufsort64.h" //用64位的，适用于全长基因组
#include "InArray.h"
#ifdef MYDEBUG
#define jj printf
#else
#define jj //printf
#endif
#define uint unsigned int
using namespace std;

// const uint N =10000;    //ref长度
// const int rLen=30, sLen=7, skip=4;
const int offrate =5;      //SA采样步长
// const int occS =128;    //Occ采样步长
const int occS =256;       //Occ采样步长
int oneSeed =0;            //1:强行指定seed个数，0:默认个数
int nrd2 =0;
const uint nrd =1000;      //read个数
const uint nFreSd =700;    //单个seed的频次限定
const int locc =2000;      //seed平均最小Occ（不含两端值）
const int hocc =3000;      //seed平均最大Occ
const uchar sm =17;        //小于等于直接跳（sm越大，Q可以越小，速度越快）
const int dp =64;          //树高（这个作用不大，太小反而使nLF增大）
// const int rLen=101, sLen=22, skip=10; //场景1：seed
const int rLen=101, sLen=15, skip=1;      //场景2：kmer
// const uint N =10000*10000*32; //ref长度
// long N =0;             //ref长度
uint nLF =0;
int nLF1=0;
// int bc =hocc;          //队列表长（最多需要hocc*nSeed）
// long dollar =0;
uchar *seq;               //含$符
uchar *BWT;               //BWT串
u32 *cBWT;                //压缩BWT串（以32为桶大小）
// long *pOcc;            //预存Occ
// long *sSA;             //SA采样
uint N1 =0;
uint dollar1 =0;
uint C1[6];                //C表
uint *pOcc1;
uint *sSA1;
struct read* readSet[nrd]; //指向read结构体的指针数组
int nSeeds;                //seed个数
// map<uint, int> *m;      //映射表
// m =new map<uint, int>;  //hash键值对：原始排位、0开头的下标
uint* m;                   //上述功能，改用数组存，提高效率 <区间起始排位、区间的大小>
uint* L;
uint* R;
// long C[6];             
int code[256];
char dna[4];
// uint maxi=0, mini=N+2;
string fseq="seq.txt.2"; //处理hg19.fa，全为大写，重复区的小写也转为大写，N和n也随机化为大写
string fbwt="bwt.txt";


string *covInfo(const string &data, int col)
{ //返回指针
    int len = data.length();
    string *ss = new string[col];
    int i = 0, num = 0;
    while (i < len)
    {
        char ch = data.at(i);
        if (ch == ' ')
            num++;
        else
            ss[num].push_back(ch);
        i++;
        if (num == col)
            break;
    }
    return ss;
}

void coding(){
    code['$']=0; code['A']=1; code['C']=2; code['G']=3; code['T']=4; //对应C表编码
    dna[0]='A'; dna[1]='C'; dna[2]='G'; dna[3]='T';
}

uint Occ0(char c, uchar BWT[], uint pos, bool fdol){ //最后一个参数:是否考虑$
    // cout<<"in..."<<"c:"<<c<<" pos:"<<pos<<endl;
    if(pos>N1+1 || pos<0){
        cout<<" (Occ0)pos is large than N, "<<pos<<endl;
        exit(1);
    }
    uint cnt[5]={0}; //存放各字符的rank值
    for(uint i=0; i<=pos; i++){ //从头到尾遍历找，是哪个字符，谁加一
        cnt[code[BWT[i]]]++; //强力扫
    }
    if(fdol && pos>dollar1) cnt[1]--; //减掉$放入的A
    // cout<<"out..."<<"code[c]:"<<code[c]<<" Occ:"<<cnt[code[c]]<<endl;
    return cnt[code[c]];
}

inline static int popCntDna(int c, u64 data){ //2bit序列的DNA_popCnt
    u64 mask[4];
    mask[0]=0xffffffffffffffff; //依次ACGT
    mask[1]=0xaaaaaaaaaaaaaaaa;
    mask[2]=0x5555555555555555;
    mask[3]=0x0000000000000000;
	u64 x0 = data ^ mask[c];
	u64 x1 = (x0 >> 1);
	u64 x2 = x1 & (0x5555555555555555);
	u64 x3 = x0 & x2;
    u64 tmp= __builtin_popcountll(x3); //内置64位popCnt
    return tmp;
}

inline void Occ(uint pos, uint cnt[], int sz, int n[], bool flag[], bool* f[4], uchar& c){ //最后一个参数:是否考虑$
    pos =pos+1;
    if(pos>N1+2 || pos<0){ cout<<" (Occ)pos is large than N, "<<pos<<endl; exit(1);}
    uint cntA=0, cntC=0, cntG=0, cntT=0;
    uint p=((pos>>8)<<2) + ((pos>>8)<<4);  //occ + bwt          
    bool tmp;
    cntA+= cBWT[p++];                     //块前
    cntC+= cBWT[p++];
    cntG+= cBWT[p++];
    cntT+= cBWT[p++];
    int lens =pos%occS;                    //块内
    u64 data=0;
    for(int i=0; i<lens>>5; i++){
        data =*(u64*)&(cBWT[p]);
        cntA+= popCntDna(0, data);         //1. 解32
        cntC+= popCntDna(1, data);
        cntG+= popCntDna(2, data);
        cntT+= popCntDna(3, data);
        p+= 2;
    }
    int lens_less32 =lens%32;
    if(lens_less32){
        data =*(u64*)&(cBWT[p]);           //取含不足32的u64块
        data =data<<((32-lens_less32)<<1);      //右移多余的位
        cntA+= popCntDna(0, data);         
        cntC+= popCntDna(1, data);
        cntG+= popCntDna(2, data);
        cntT+= popCntDna(3, data);
        cntA =cntA-(32-lens_less32);       //去掉移位产生的00
    }
    if(pos>(dollar1+1)) cntA--;             //减掉$放入的A
    cnt[1]=cntA; cnt[2]=cntC; cnt[3]=cntG; cnt[4]=cntT;
    c =(cBWT[(((pos-1)>>8)<<2) + 4 + ((pos-1)>>4)] >> (((pos-1)%16)<<1)) &0x3;

    if(sz){                                //需要继续向后，返回f、n
        n[0]=1; n[1]=1; n[2]=1; n[3]=1;    //n是各分组个数，cnt为各rank值
        for(int i=0; i<4; ++i){            //生成f，并将首位置1
            // f[i] =new bool[sz+1];          //f中首位放是否含0，n[]比实际大1
            f[i][0] =1;                    //标记位，代表为全1
        }
        int k=0, j=1;                          //变量复用，减少开销
        for(uint i=pos; i<pos+sz; ++i){    //count, locate的sz为0，不会执行这里
            // k= (cBWT[i>>4]<<((15-i%16)*2))>>30; //先取到值
            tmp =flag[j++];
            // k= (cBWT[i>>4]>>((i%16)<<1))&0x3;
            k= (cBWT[((i>>8)<<2) + 4 + (i>>4)] >> ((i%16)<<1)) &0x3;      //bwt串上取值
            f[k][0] =f[k][0] &&tmp;     //一旦有元素为0，则将f的首位置0（首位有0则为0，初值为1）
            f[k][n[k]++] =tmp;          //产生各组的flag
        }
    }
}

inline void Occ(uint pos, uint cnt[], int sz, int n[], bool flag[], bool* f[4], uint o[], uint* oo[4]){ //最后一个参数:是否考虑$
    pos =pos+1;
    if(pos>N1+2 || pos<0){ cout<<" (Occ2)pos is large than N, "<<pos<<endl; exit(1);}
    uint cntA=0, cntC=0, cntG=0, cntT=0;
    uint p=((pos>>8)<<2) + ((pos>>8)<<4);             //occ + bwt 
    cntA+= cBWT[p++];                                 //块前
    cntC+= cBWT[p++];
    cntG+= cBWT[p++];
    cntT+= cBWT[p++];
    int lens =pos%occS;                               //块内
    u64 data=0;
    for(int i=0; i<lens>>5; i++){
        data =*(u64*)&(cBWT[p]);
        cntA+= popCntDna(0, data);                    //1. 解32
        cntC+= popCntDna(1, data);
        cntG+= popCntDna(2, data);
        cntT+= popCntDna(3, data);
        p+= 2;
    }
    int lens_less32 =lens%32;
    if(lens_less32){
        data =*(u64*)&(cBWT[p]);                      //取含不足32的u64块
        data =data<<((32-lens_less32)<<1);            //右移多余的位
        cntA+= popCntDna(0, data);         
        cntC+= popCntDna(1, data);
        cntG+= popCntDna(2, data);
        cntT+= popCntDna(3, data);
        cntA =cntA-(32-lens_less32);                  //去掉移位产生的00
    }
    if(pos>(dollar1+1)) cntA--;                       //减掉$放入的A
    cnt[1]=cntA; cnt[2]=cntC; cnt[3]=cntG; cnt[4]=cntT;
    if(sz){                                           //需要继续向后，返回f、n
        n[0]=1; n[1]=1; n[2]=1; n[3]=1;               //n是各分组个数，cnt为各rank值
        f[0][0]=1; f[1][0]=1; f[2][0]=1; f[3][0]=1;   //标记位，代表为全1
        int k=0, j=1; bool tmp;                       //变量复用，减少开销
        for(uint i=pos; i<pos+sz; ++i){               //count, locate的sz为0，不会执行这里
            tmp =flag[j++];                           //j为pos后bwt串上的移动指针
            k= (cBWT[((i>>8)<<2) + 4 + (i>>4)] >> ((i%16)<<1)) &0x3;      //bwt串上取值
            f[k][0] =f[k][0] &&tmp;                   //一旦有元素为0，则将f的首位置0（首位有0则为0，初值为1）
            f[k][n[k]++] =tmp;                        //产生各组的flag（同时计数++）
            oo[k][n[k]-2] =o[j-2];                    //传递各组中的o
        }
    }
}

void save(string fname, uchar array[], uint n){    //存seq和bwt用
    ofstream of;
    of.open(fname.c_str());
    long cnt[256]={0};                             //统计字符数
    int nChar =0;
    for(uint i=0; i<n; i++){
        of<<array[i];                              //逐字符的写
        cnt[array[i]]++;                           //计数
    }
    for(int i=0; i<256; i++){
        if(cnt[i] >0) {
            nChar++;
            cout<<i<<" "<<(char)i<<": "<<cnt[i]<<endl;
        }
    }
    cout<<"charNum:"<<nChar<<endl;
    of.close();
}

struct group
{
    uint l, r;
    int step; bool* flag;
    uint* o;                                                       //glocate+用
    group(){}                                                      //默认构造函数
    inline group(uint l, uint r, int step, bool flag[]){           //glocate用（不带o数组）
        this->l=l;                                                 //必须要有this
        this->r=r;
        this->step=step;
        this->flag=new bool[r-l+2];                                //首个为是否含0标记
        memcpy(this->flag, flag, r-l+2);
    }
    inline group(uint l, uint r, int step, bool flag[], uint o[]){ //glocate+用
        this->l=l;                                                 //必须要有this
        this->r=r;
        this->step=step;
        this->flag=new bool[r-l+2];                                //首个为是否含0标记
        this->o=new uint[r-l+1];                                   //new分配空间，把传入的o放入new的空间中
        // for(int i=0; i<r-l+1; i++){
        //     this->flag[i]=flag[i];
        //     this->o[i]=o[i];
        // }
        // this->flag[r-l+1]=flag[r-l+1];                          //处理多出来的一个
        memcpy(this->flag, flag, r-l+2);
        memcpy(this->o, o, (r-l+1)*4);
    }
};

class MyQueue{                           //循环队列
    int f, r, c;                         //首尾下标
    group *array;                        //数组

    public:
    MyQueue(int capacity){               //构造函数
        this->c =capacity;
        this->f =0;
        this->r =0;
        this->array =new group[capacity]; //调用默认构造函数
    }
    inline void enQueue(group e){         //入队
        array[r] =e;
        r =(r+1)%c;
        if(r ==f){
            cout<<"Error: queue is full. Exit..."<<endl;
            exit(1);
        }
    }
    inline group* deQueue(){              //出队
        group* e =&array[f];
        f =(f+1)%c;
        return e;
    }
    inline int size(){
        return (r-f+c)%c;                 //(尾-头+表长)%表长
    }
};

inline bool nALL1(bool flag[], int& n){
    // int len=sizeof(flag)/sizeof(int);
    // cout<<"nALL1...len:"<<len<<endl;
    for(int i=0; i<n; ++i){
        // cout<<flag[i];
        if(!flag[i]){
            return true; //非全1
        }
    }
    // cout<<endl;
    return false; //全1
}

// void genRef0(){
//     char* ss[5]={"AGACGT","GGACGT","AGACGT","GGACGT","AGACGT"}; //含patter的字串
//     int tmpC=0; //构建测试序列
//     for(int i=0; i<N; i++){
//         tmpC=lrand48()&0x3; //产生0123
//         switch (tmpC)
//         {
//             case 0:
//                 ref[i]='A';
//                 break;
//             case 1:
//                 ref[i]='C';
//                 break;
//             case 2:
//                 ref[i]='G';
//                 break;
//             case 3:
//                 ref[i]='T';
//                 break;
//             default:
//                 break;
//         }
//         // cout<<ref[i];
//     }
//     int tmp_pos, b=990, a=10;
//     for(int i=0; i<5; i++){
//         tmp_pos =(rand()%(b-a+1))+ a; //产生5个随机位置，10~990之间
//         // cout<<tmp_pos<<" ";
//         // cout<<strlen(ss[i])<<" ";
//         // cout<<ss[i]<<endl;
//         for(int j=0; j<strlen(ss[i]); j++){
//             ref[tmp_pos+j] =ss[i][j]; //用ss序列，替换 
//         }
//     }
// }

void genRef(){ //人工合成数据
    seq =new uchar[N1+1];
    // char* ss[5]={"AGACGT","GGACGT","AGACGT","GGACGT","AGACGT"};
    const char *ss[5]={"AGACGTAGACGTAGACGT","GGACGTGGACGTGGACGT","AGACGTAGACGTAGACGT","GGACGTGGACGTGGACGT","AGACGTAGACGTAGACGT"}; //含patter的字串（指向常量的指针
    int tmpC=0; //构建测试序列2（邻近重复
    for(uint i=0; i<N1; i++){
        tmpC=lrand48()&0x3; //产生0123
        switch (tmpC)
        {
            case 0:
                seq[i]='A';
                break;
            case 1:
                seq[i]='C';
                break;
            case 2:
                seq[i]='G';
                break;
            case 3:
                seq[i]='T';
                break;
            default:
                break;
        }
        // cout<<ref[i];
    }
    uint tmp_pos, b=N1-20, a=20; //替换的范围
    for(int i=0; i<5; i++){
        tmp_pos =(rand()%(b-a+1))+ a; //产生5个随机位置，20~980之间
        // cout<<tmp_pos<<" ";
        // cout<<strlen(ss[i])<<" ";
        // cout<<ss[i]<<endl;
        for(int j=0; j<strlen(ss[i]); j++){
            seq[tmp_pos+j] =ss[i][j]; //用ss序列，替换 
        }
    }
    seq[N1] ='$';
    // seq[N] ='A'; //为了便于2bit编码
    save(fseq, seq, N1+1); //存seq
}

uint readHg19(string file){ //读入hg19参考基因组, 生成seq.txt序列
    int myseed = 11;
	srand48(myseed); //随机种子
    seq =new uchar[4000000000]; //预开辟40亿长度
    // string file ="/home/lab/gll/data_aln/bowtie2-master/example-hg19/reference/hg19.fa";
    // string file ="/home/lab/gll/data_aln/bowtie2-master/example/reference/lambda_virus.fa"; //测试
    ifstream in(file.c_str());
    if(in.fail()){cout<<"file open error!"<<endl; exit(1);}
    string data;
    uint k=0;
    while(getline(in, data)){ //整行读
        // cout<<data<<endl;
        if(data.length()==0 || data.at(0)=='>') continue; //跳过最后空行，>染色体ID行
        for(int i=0; i<data.length(); i++){
            char c =data.at(i);
            if(c=='N' || c=='n'){
                c=dna[lrand48()&3]; //将N随机化
            }
            else if(c>'Z'){ //重复区的小写字符
                c =c-32; //转为大写
            }
            seq[k++]=c; //放入序列中
        }
    }
    // seq[k]='$';
    // // save(fseq, seq, k+1); //存seq
    
    // save(fseq, seq, k); //存seq (不含$, 2024.3.14测FMI建索引用) seq.txt
    return k; //实际长度
}

// int genRead_Seeds(){
    // int nSeeds=(rLen-sLen)/skip +1;
    // seeds=new uchar* [nSeeds]; //二维数组分配内存
    // for(int i=0; i<nSeeds; i++){
    //     seeds[i]=new uchar[sLen];
    // }
    // uint pos, b=N-120, a=0;
    // pos=(rand()%(b-a+1))+ a;
    // cout<<"reads pos:"<<pos<<" n:"<<nSeeds<<endl;
    // for(int i=0; i<nSeeds; i++){
    //     for(int j=0; j<sLen; j++){
    //         seeds[i][j] =seq[pos+i*skip+j]; //从seq上读，如果直接载入索引，seq无值报错
    //     }
    // }
    // for(int i=0; i<nSeeds; i++){
    //     for(int j=0; j<sLen; j++){
    //         cout<<seeds[i][j];
    //     }
    //     cout<<endl;
    // }
    // delete[] seq;
    // return nSeeds;
// }

struct read{ //保存read信息
    bool rep; //重复区
    uint pos, *L, *R;
    char **seeds;
    read(int nSeeds){
        seeds =new char* [nSeeds]; //二维数组分配内存
        for(int i=0; i<nSeeds; i++){
            seeds[i]=new char[sLen];
        }
        L =new uint[nSeeds]; //存放SA区间
        R =new uint[nSeeds];
    }
    void delRead(int nSeeds){ //需要有析构函数
        delete[] L;
        delete[] R;
        for(int i=0; i<nSeeds; i++){
            delete[] seeds[i];
        }
        delete[] seeds;
    }
};

void count(char seed[], int sLen, uint& L, uint& R, bool prt){
    int nLF =0; uchar unUse;//统计LF次数
    char c=seed[sLen-1];
    L=C1[code[c]]; //首字符查C表
    R=C1[code[c] +1];
    for(int i=sLen-2; i>=0; i--){ 
        c=seed[i];
        uint cnt[5]; int n[4]; bool flag[1]; bool* f[4];
        // cout<<"bef... i:"<<i<<" c:"<<c<<" L:"<<L<<" C:"<<C[code[c]]<<" code:"<<code[c]<<endl;
        // for(int i=0; i<6; i++){cout<<C[i]<<" ";} cout<<endl; //打印C表
        // Occ(cBWT, L-1, cnt, 0, n, flag, f, unUse);
        Occ(L-1, cnt, 0, n, flag, f, unUse);
        L=C1[code[c]] +cnt[code[c]];
        // cout<<"aft... i:"<<i<<" c:"<<c<<" L:"<<L<<endl;

        // Occ(cBWT, R, cnt, 0, n, flag, f, unUse);
        Occ(R, cnt, 0, n, flag, f, unUse);
        R=C1[code[c]] +cnt[code[c]] -1;
        nLF++;
    }
    if(prt){
        // cout<<"count (nLF:"<<nLF<<endl;
        cout<<"seed:"; for(int i=0; i<sLen; i++){cout<<seed[i];}
        cout<<" count:"<<R-L+1<<" L:"<<L<<" R:"<<R<<endl;
    }
}

int genRead_Seeds_Mmap(uint nrd){
    char *seq; struct stat stat;
    int fd =open("hg19.seq", O_RDWR, 0);
    if(fd <0) {cout<<"file open error!"<<endl;}
    fstat(fd, &stat);
    seq =(char*)mmap(NULL, stat.st_size, PROT_READ, MAP_SHARED, fd, 0); //mmap使用

    int nSeeds=(rLen-sLen)/skip + ((((rLen-sLen)%skip)==0)?0:1); //有余数则加一，无余数则加零
    if(oneSeed==1){
        nSeeds=1; //为调试方便，强行指定seed为1个
    }

    uint pos, b=N1-120, a=0;

    // for(int k=0; k<nrd; k++){
    int k=0; bool stop;
    while(true){
        bool sdCheck =true; //检查单个seed是否满足nFreSd阈值
        stop =0;
        if(k==nrd) break; //找到足够的read，则跳出循环
        bool rep=false;
        pos=(rand()%(b-a+1))+ a; //随机选取的read位置
        struct read* tmp= new struct read(nSeeds);
        for(int i=0; i<nSeeds; i++){
            for(int j=0; j<sLen; j++){
                char c =seq[pos+i*skip+j];
                if(c=='N'||c=='n'){ //碰到N或n，丢掉这个read
                    stop =1;
                }
                if(c>'Z'){
                    c =c-32;         //重复区的小写，转为大写    
                    rep =true;
                }
                tmp->seeds[i][j] =c;     //从seq上读
            }
            for(int j=0; j<i; j++){      //检查是否存在相同pattern
                if(strcmp(tmp->seeds[i],tmp->seeds[j])==0){   
                    stop =1;
                    break;
                }     
            }
            if(stop){break;}
        }
        if(stop){continue;}
        tmp->rep =(rep)?true:false;
        tmp->pos =pos;
        // if(tmp->rep){
        //     cout<<"reads# "<<k<<"  pos:"<<pos<<" n:"<<nSeeds<<endl;
        // }else{
        //     cout<<"reads. "<<k<<"  pos:"<<pos<<" n:"<<nSeeds<<endl;
        // }
        int sz=0;
        // cout<<"k:"<<k<<endl;
        for(int j=0; j<nSeeds; j++){
            // cout<<"pos:"<<pos; cout<<" seed:"; for(int i=0; i<sLen; i++){cout<<tmp->seeds[j][i];}
            count(tmp->seeds[j], sLen, tmp->L[j], tmp->R[j], false);                                                     //count查询
            // cout<<" j:"<<j<<" l:"<<tmp->L[j]<<" r:"<<tmp->R[j]<<" sz:"<<tmp->R[j]-tmp->L[j]+1<<endl;
            if(tmp->R[j]-tmp->L[j]+1 <nFreSd){                 //限制单个pattern频次
                sdCheck =false;
                break;
            }
            sz+= tmp->R[j]- tmp->L[j] +1;
        }
        if(sdCheck && sz> nSeeds*locc && sz< nSeeds*hocc){                                                               //出现次数阈值（不包含两端值）
            readSet[k++] =tmp;                                                                                           //满足条件的read，则k++
            if(tmp->rep){
                cout<<"reads# "<<k<<"  pos:"<<setw(12)<<pos<<"  n:"<<nSeeds<<"  sz:"<<sz<<"  ave:"<<sz/nSeeds/1.0<<endl;
            }else{
                cout<<"reads. "<<k<<"  pos:"<<setw(12)<<pos<<"  n:"<<nSeeds<<"  sz:"<<sz<<"  ave:"<<sz/nSeeds/1.0<<endl;
            }
        }else{
            tmp->delRead(nSeeds);                                                                                        //用不到的reads清除
            delete tmp;
        }
    }

    munmap(seq, stat.st_size);
    return nSeeds;
}

void inData(u32 *data, uint pos, int v){
    data[pos/16] =data[pos/16]|(v<<((pos%16)*2));
}

void genC_BWT_sSA(){
    // uchar mask[4][4] ={{0xFC,0xF3,0xCF,0x3F}, //[字符][位置]
    //                    {0xFD,0xF7,0xDF,0x7F},
    //                    {0xFE,0xFB,0xEF,0xBF},
    //                    {0x,0x,0x,0x}};
    // cBWT[i/4] =cBWT[i/4] | mask[code[seq[SA[i] -1]]-1][i%4];
    // cBWT =new InArray((fm_int)N+1, 2); //2bit编码
    cBWT =new u32[(N1+1)/16 + (N1+1)/occS*4 +1]; //2bit编码
    memset(cBWT,0,((N1+1)/16 + (N1+1)/occS*4 +1)*4);
    BWT =new uchar[N1+1];
    // pOcc1 =new uint[((N1+1)/occS +1)*4]; //预存Occ
    long *SA =new long[N1+1]; //SA数组
    divsufsort64(seq, SA, N1+1); //产生SA数组

    uint nA=0, nC=0, nG=0, nT=0, k=0, nOcc=0;
    for (uint i=0; i<N1+1; i++) //产生BWT串, 统计字符频率
    {
        if(i%occS ==0){ //pOcc预存
            uint tmp=(i/16)+(nOcc*4);
            cBWT[tmp++]=nA;
            cBWT[tmp++]=nC;
            cBWT[tmp++]=nG;
            cBWT[tmp]=nT;
            nOcc++;
        }
        if (SA[i] == 0){
            dollar1=i;
            cout<<"dollar:"<<dollar1<<endl;
            BWT[i]='A'; //'$' 36
            // cBWT->SetValue((fm_int)i, (u64)code['G']-1);
            // inData(cBWT,i,code['A']-1);
            inData(cBWT, i+(nOcc*64), code['A']-1);
            // cout<<"sa0:"<<i<<" $:"<<ref[N]<<" BWT[]:"<<BWT[i]<<endl;
        }else{
            BWT[i] =seq[SA[i] -1];
            inData(cBWT, i+(nOcc*64), code[seq[SA[i]-1]] -1); //先存pOcc，再放bwt
        }
        switch (BWT[i]) //C表统计
        {
            case 'A':
                nA++;
                break;
            case 'C':
                nC++;
                break;
            case 'G':
                nG++;
                break;
            case 'T':
                nT++;
                break;
            default:
                break;
        }
    }
    // save(fbwt, BWT, N+1); //存bwt
    
    nA--; //减掉$放入的A
    C1[0]=0; C1[1]=1; C1[2]=nA+C1[1]; C1[3]=nC+C1[2]; C1[4]=nG+C1[3]; C1[5]=nT+C1[4];
    // C[0]=0; C[1]=0; C[2]=nA+C[1]; C[3]=nC+C[2]; C[4]=nG+C[3]; C[5]=nT+C[4]; //不含$
    for(int i=0; i<6; i++){
        cout<<C1[i]<<" ";
    }
    // getchar();
    sSA1 =new uint[((N1+1)>>offrate) +1]; //>>和+，一定注意加括号，+的优先级更高
    for(uint i=0,j=0; i<N1+1; i++){
        if(i%(1<<offrate)==0){ //下标采样
            sSA1[j]=SA[i];
            j++;
        }
    }
    if(true){ //索引存磁盘
        string file ="glocate.idx.hg19.32.256.cross";
        fstream fs(file.c_str(),ios::out|ios::binary);
        if(!fs.is_open()){cout<<"file is not exits!"<<endl;}
        fs.write((char*)&N1, sizeof(uint)); //序列长度（不含$
        fs.write((char*)&dollar1,sizeof(uint)); //dollar位置
        fs.write((char*)C1, sizeof(uint)*6); //C表
        // fs.write((char*)pOcc1, sizeof(uint)*(((N1+1)/occS +1)*4)); //Occ预存
        fs.write((char*)cBWT, sizeof(u32)*((N1+1)/16 + (N1+1)/occS*4 +1)); //2bit的BWT串
        fs.write((char*)sSA1, sizeof(uint)*(((N1+1)>>offrate) +1)); //采样SA（这里>>和+，要注意加上括号）
        fs.close();
    }
    delete[] SA; //后续计算基于sSA
    // delete[] seq; //放到genRead_Seeds中清
    delete[] BWT; //后续计算基于cBWT
}

void loadIdx(string file){ //载入索引
    // string file ="glocate.idx.hg19.32.256.cross";
    fstream fs(file.c_str(),ios::in|ios::binary);
    if(!fs.is_open()){cout<<"file is not exits!"<<endl;}
    fs.read((char*)&N1, sizeof(uint)); //序列长度（不含$
    cout<<"N1:"<<N1<<endl;
    fs.read((char*)&dollar1, sizeof(uint)); //dollar位置
    cout<<"dollar1:"<<dollar1<<endl;
    fs.read((char*)C1, sizeof(uint)*6); //C表
    cout<<"C1:"; for(int i=0; i<6; i++){cout<<C1[i]<<" ";} cout<<endl;
    // pOcc1 =new uint[((N1+1)/occS +1)*4];
    // fs.read((char*)pOcc1, sizeof(uint)*(((N1+1)/occS +1)*4)); //Occ预存
    cBWT =new u32[(N1+1)/16 + (N1+1)/occS*4 +1];
    fs.read((char*)cBWT, sizeof(u32)*((N1+1)/16 + (N1+1)/occS*4 +1)); //2bit的BWT串

    // cout<<"cBWT finish!"<<endl;
    // for(uint i=0; i<(N1+1)/16 + (N1+1)/occS*4 +1; i++){ //查看cBWT的内容
    //     if(i%20==0) getchar();
    //     cout<<cBWT[i]<<" ";
    // }

    sSA1 =new uint[(N1+1)/(1<<offrate) +1];
    fs.read((char*)sSA1, sizeof(uint)*((N1+1)/(1<<offrate) +1)); //采样SA

    cout<<"load finish!"<<endl;
}

// void loadIdx(){ //载入索引
//     string file ="glocate.idx.hg19.32";
//     fstream fs(file.c_str(),ios::in|ios::binary);
//     if(!fs.is_open()){cout<<"file is not exits!"<<endl;}
//     fs.read((char*)&N, sizeof(long)); //序列长度（不含$
//     cout<<"N:"<<N<<endl;
//     fs.read((char*)&dollar,sizeof(long)); //dollar位置
//     cout<<"dollar:"<<dollar<<endl;
//     fs.read((char*)C, sizeof(long)*6); //C表
//     cout<<"C:"; for(int i=0; i<6; i++){cout<<C[i]<<" ";} cout<<endl;
//     pOcc =new long[((N+1)/occS +1)*4];
//     fs.read((char*)pOcc, sizeof(long)*(((N+1)/occS +1)*4)); //Occ预存
//     cBWT =new u32[(N+1)/16 +1];
//     fs.read((char*)cBWT, sizeof(u32)*((N+1)/16 +1)); //2bit的BWT串
//     sSA =new long[(N+1)/(1<<offrate) +1];
//     fs.read((char*)sSA, sizeof(long)*((N+1)/(1<<offrate) +1)); //采样SA

//     cout<<endl;
//     N1=N;           cout<<"N1:"<<N1<<endl;
//     dollar1=dollar; cout<<"dollar1:"<<dollar1<<endl;
//     cout<<"C1:";    for(int i=0; i<6; i++){C1[i]=C[i]; cout<<" "<<C1[i];} cout<<endl;
//     pOcc1 =new uint[((N1+1)/occS +1)*4]; for(uint i=0; i<((N1+1)/occS +1)*4; i++){pOcc1[i]=pOcc[i];}
//     sSA1 =new uint[(N1+1)/(1<<offrate) +1]; for(uint i=0; i<(N1+1)/(1<<offrate) +1; i++){sSA1[i]=sSA[i];}
//     delete[] pOcc;
//     delete[] sSA;
//     cout<<"load finish!"<<endl;
// }

uint* locate(uint L[], uint R[], int n, int& sz){
    sort(&L[0],&L[n]); //排序
    sort(&R[0],&R[n]); //排序
    uint cnt[5]; int nn[4]; bool flag[1]; bool* f[4];
    nLF=0; nLF1=0;
    sz=0; uchar c;
    for(int i=0; i<n; i++){
        sz +=(R[i]-L[i]+1);
    }
    for(int i=0; i<4; ++i){                                          //避免在occ中频繁new
        f[i] =new bool[sz+1];
    }
    uint* loc =new uint[sz];
    memset(loc,-1,sizeof(loc));
    int num=0;
    for(int s=0; s<n; s++){
        for(uint k=L[s]; k<=R[s]; k++){
            int step=0;
            uint i=k;
            while(i%(1<<offrate) !=0){
                Occ(i, cnt, 0, nn ,flag, f, c);
                step++; nLF++; //LF计数
                i=C1[c+1] +cnt[c+1] -1;
            }
            loc[num++]= sSA1[i>>offrate] +step;
        }
    }
    // cout<<"locate  (nLF:"<<nLF<<endl;                             //以下为显示locate结果
    // sort(&loc[0],&loc[sz]);                                       //排序
    // for(int i=0; i<sz; i++){cout<<loc[i]<<" ";} cout<<endl;
    return loc;                                                      //返回位置个数
}

// void bwtcpy(u32 *ctbwt, uchar tbwt[], u32 *cBWT, int l, int r){
void bwtcpy(uchar tbwt[], u32 *cBWT, uint l, uint r){
    for(uint i=l,j=0; i<=r; i++,j++){
        int v= (cBWT[i/16]<<((15-i%16)*2))>>30; //先取到值
        // inData(ctbwt, j, v); //放入值
        tbwt[j] =dna[v];
    }
}

uint* glocate(uint L[], uint R[], int n){
    sort(&L[0],&L[n]); //排序
    sort(&R[0],&R[n]); //排序
    nLF=0; nLF1=0;
    bool *flag; bool* f[4];
    int sz=0, num=0, rt, nn[4], nn2[4], step=0, tmp=0;
    uint l=0, r=0, ii, cnt[5], cnt2[5];
    uchar c=0;
    struct group* pg;
    MyQueue Q(hocc*n);                                              //队列大小，默认60000，每个seed分配的表长，再乘seed个数
    for(int i=0; i<n; ++i){
        l=L[i]; r=R[i];
        sz +=r-l+1;
        bool zero[r-l+2]={0};                                    //首个为是否含0
        Q.enQueue(group(l,r,0,zero));                            //首次入队
    }
    for(int i=0; i<4; ++i){                                      //避免在occ中频繁new
        f[i] =new bool[sz+1];
    }
    rt =1<<offrate;
    uint* loc =new uint[sz];                                     //栈空间
    memset(loc,-1,sizeof(loc));
    while(num <sz){                                              //执行R-L+1次
        pg =Q.deQueue();                                         //取队头元素
        l=pg->l; r=pg->r; step=pg->step; flag=pg->flag;
        for(uint i=l; i<=r; ++i){
            if(flag[i-l+1] ==1) continue;                        //flag==1，当前i是已解的
            if(i%rt ==0){                                        //当前i未解，且为采样点
                loc[num++] =sSA1[i>>offrate]+step;
                flag[i-l+1]=1;                                   //true表已解
            }
        }
        Occ(l-1, cnt, r-l+1, nn, flag, f, c);                    //bwt区间字符计数，flag更新，都在Occ内完成
        nLF++; step++; nLF1++;
        for(int k=0; k<4; ++k){
            if( nn[k]>1 && (f[k][0]==0) ){                       //含未解的，需要继续找，才会加入
                l=C1[k+1]+cnt[k+1];
                r=l+nn[k]-2;
                if( ((nn[k]-1)<=sm) || (step>dp) ){              //元素少直接找
                    for(uint i=l; i<=r; ++i){
                        if(f[k][i-l+1]) continue;                //已解
                        ii =i; 
                        tmp =step;
                        while(ii%rt !=0){
                            Occ(ii, cnt2, 0, nn2, flag, f, c);
                            nLF++; tmp++; ii =C1[c+1] +cnt2[c+1] -1;
                        }
                        loc[num++] =sSA1[ii>>offrate] +tmp;
                    }
                }else{
                    Q.enQueue(group(l, r, step, f[k]));
                }
            }
        }
    }
    // cout<<"glocate (nLF:"<<nLF<<endl; //以下为显示glocate结果
    // sort(&loc[0],&loc[sz]); //排序
    // for(int i=0; i<sz; i++){cout<<loc[i]<<" ";} cout<<endl;
    return loc;
}

// inline int mp(uint& o, uint m[], int& n){ //给排位，得0起始的下标.（见2023.9.1周报，三种方案中的线性扫）
//     uint t=0;
//     int k=0, c=0, s=0;
//     for(int i=0; i<n; i++){   //对n个区间扫描
//         t=o-m[i<<1];          //o-开始排位
//         s=m[(i<<1) +1];       //取这个区间的个数
//         if(t>=0 && t<s){      //判断在这个区间中
//             k=c+t;            //之前的，加当前所处区间
//             break;            //找到则退出扫描
//         }
//         c=c+s;                //累加之前的个数
//     }
//     // cout<<"o:"<<k<<endl;
//     // getchar();
//     return k;
// }

// inline static bool in(uint& o){ //判断i排位，是否在这几个非连续区间段内
//     for(int i=0; i<nSeeds; ++i){
//         if(o>=m[i<<1] && (o<(m[i<<1]+m[(i<<1)+1]))){
//             return true;
//         }
//     }
//     return false;
// }

// inline static int mp(uint& o){ //给排位(一定在区间内)，得0起始的下标.（见2023.9.1周报，三种方案中的线性扫）
//     int k=0;
//     for(int i=0; i<nSeeds; ++i){   //对n个区间扫描
//         if(o>=m[i<<1] && (o<(m[i<<1]+m[(i<<1)+1]))){      //判断在这个区间中
//             k=k+o-m[i<<1];            //之前的，加当前所处区间
//             return k;            //找到则退出扫描
//         }
//         k=k+m[(i<<1)+1];                //累加之前的个数
//     }
//     // cout<<"o:"<<k<<endl;
//     // getchar();
//     return k;
// }

inline bool in(uint& o){ //判断i排位，是否在这几个非连续区间段内
    for(int i=0; i<nSeeds; ++i){
        if(o>=L[i] && o<=R[i]){
            return true;
        }
    }
    return false;
}

// inline static int mp(uint& o){ //给排位(一定在区间内)，得0起始的下标.（见2023.9.1周报，三种方案中的线性扫）
//     // if(o>maxi || o<mini){ return -1;}
//     for(int i=0; i<nSeeds; ++i){                          //对n个区间扫描
//         if(o>=L[i] && o<=R[i]){      //判断在这个区间中
//             return o-L[i]+m[i];                      //之前的，加当前所处区间
//         }
//     }
//     return -1;
// }

inline static int mp(uint& o){       //折半查
    int l=0, h=nSeeds-1, mid;
    while(l<=h){
        mid=(l+h)/2;
        if(o>=L[mid] && o<=R[mid]){
            return o-L[mid]+m[mid];
        }else if(o>=L[mid]){         //在右
            l=mid+1;
        }else{                       //在左
            h=mid-1;
        }
        // cout<<"mid:"<<mid<<" l:"<<l<<" h:"<<h<<endl;  getchar();
    }
    return -1;
}

uint* glocatePlus(uint L[], uint R[], int n){
    sort(&L[0],&L[n]); //排序
    sort(&R[0],&R[n]); //排序
    nLF=0; nLF1=0;
    bool *flag; bool* f[4];
    uint *o; uint* oo[4];  
    int sz=0, num=0, rt, nn[4], nn2[4], step=0, tmp=0, t=0, h=0, k=0;
    uint l=0, r=0, ii, cnt[5], cnt2[5];
    uchar c=0;
    struct group* pg;                                                    //指针
    MyQueue Q(hocc*n);                                                   //循环队列
    for(int i=0; i<n; ++i){                                              //依次处理n个seeds
        m[i]=sz;                                                         
        sz +=(R[i]-L[i]+1);
        bool zero[R[i]-L[i]+2]={0};
        uint o[R[i]-L[i]+1]={0};                                         //又定义局部量，前面有uint*
        for(uint j=0; j<R[i]-L[i]+1; ++j){
            o[j]=L[i]+j;                                                 //产生初始状态的o数组
        }                  
        Q.enQueue(group(L[i],R[i],0,zero,o));                            //首次入队（比glocate多o数组）
    }
    for(int i=0; i<4; ++i){                                              //避免在occ中频繁new
        f[i] =new bool[sz+1];
        oo[i] =new uint[sz];
    }
    rt =1<<offrate;
    uint* loc =new uint[sz];                                             //已解位置
    int bef[sz];                                                         //跳了几步
    int who[sz];                                                         //subs是跳到的排位，value是当前待解排位
    memset(loc,0,sizeof(uint)*sz);
    memset(bef,-1,sizeof(int)*sz);
    memset(who,-1,sizeof(int)*sz);
    while(num <sz){                                                      //保证每个都解
        pg =Q.deQueue();                                                 //取队头元素
        l=pg->l; r=pg->r; o=pg->o; step=pg->step; flag=pg->flag;         //取出该组的信息
        for(uint i=l; i<=r; ++i){                                        //检查区间是否含采样
            if(flag[i-l+1] ==1) continue;                                //flag==1，当前i是已解的
            if(i%rt ==0){
                loc[mp(o[i-l])] =sSA1[i>>offrate]+step;
                flag[i-l+1]=1;                                           //已解（真）
                num++;
            }else if(step>0){                                            //第一步肯定在区间中，所以要加条件 (Todo:整个区间判断, 15,1,10,100时出现)
                h =mp(i);
                if(h<0) continue;                                        //检查是否跳回到原始区间
                t =mp(o[i-l]);                                           //i在所有区间中，是第几个
                bef[t] =step;                                            //bef待解排位
                who[h] =t;                                               //who碰到的排位
                flag[i-l+1]=1;                                           //已解（假）
                num++;
            }
        }
        Occ(l-1, cnt, r-l+1, nn, flag, f, o, oo);                        //不涉及o数组（flag为传入，f为传出）
        nLF++; step++; nLF1++;                                           //tmp置0
        for(int k=0; k<4; ++k){
            if( nn[k]>1 && (f[k][0]==0) ){                               //该分组含未解  
                l =C1[k+1]+cnt[k+1]; 
                r =l+nn[k]-2;
                if( ((nn[k]-1)<=sm) || (step>dp) ){                      //元素少直接找
                    for(uint i=l; i<=r; ++i){
                        if(f[k][i-l+1]) continue;                        //该元素已解，直接跳过
                        ii =i; 
                        tmp =step;
                        t =mp(oo[k][i-l]);
                        while(ii%rt !=0){
                            h =mp(ii);
                            if(h>=0){
                                bef[t]=tmp;                              //bef为待解排位的下标
                                who[h]=t;                                //who为碰到排位的下标
                                break;
                            }
                            Occ(ii, cnt2, 0, nn2, flag, f, c);           //注意，这里的cnt和nn一不小心就把前面的cnt覆盖
                            nLF++; tmp++; ii =C1[c+1]+cnt2[c+1]-1;       //更新ii
                        }
                        if(ii%rt ==0){
                            loc[t] =sSA1[ii>>offrate] +tmp;
                        }
                        num++;
                    }
                }else{
                    Q.enQueue(group(l, r, step, f[k], oo[k]));                                      
                }    
            }
        }
    }
    t=h=-1;
    for(int j=0; j<sz; ++j){
        if(loc[j] !=0){               //表明j已解
            h=j;
            while(who[h] !=-1){       //该下标可以解其它
                t=who[h];             //h能去解谁，即t
                loc[t]=loc[h]+bef[t]; //t是待解的，h是已解的下标
                who[h]=-1;            //已用，使失效
                h=t;                  //传递下去，新解的t能去解谁
            }
        }
    }
    // cout<<"gloca++ (nLF:"<<nLF<<endl;
    // sort(&loc[0],&loc[sz]); //排序
    // for(int i=0; i<sz; i++){cout<<loc[i]<<" ";} cout<<endl;
    // getchar();
    return loc;
}

void testOcc(){
    for(uint i=0; i<N1+1; i++){
        uint cnt[5]; int n[5]; bool flag[1]; bool* f[4]; uchar unUse;
        // Occ(cBWT, i, cnt, 0, n, flag, f, unUse);
        Occ(i, cnt, 0, n, flag, f, unUse);
        // for(int i=0; i<4; i++){
        //     delete[] f[i];
        // }

        char c='A';
        if(Occ0(c, BWT, i, true)!=cnt[1]){ //'A'
            cout<<"i:"<<i<<" c:"<<c<<" unEqual!"<<" Occ0:"<<Occ0(c, BWT, i, true)<<" Occ:"<<cnt[1]<<endl;
            getchar();
        }
        c='C';
        if(Occ0(c, BWT, i, true)!=cnt[2]){ //'C'
            cout<<"i:"<<i<<" c:"<<c<<" unEqual!"<<" Occ0:"<<Occ0(c, BWT, i, true)<<" Occ:"<<cnt[2]<<endl;
            getchar();
        }
        c='G';
        if(Occ0(c, BWT, i, true)!=cnt[3]){ //'G'
            cout<<"i:"<<i<<" c:"<<c<<" unEqual!"<<" Occ0:"<<Occ0(c, BWT, i, true)<<" Occ:"<<cnt[3]<<endl;
            getchar();
        }
        c='T';
        if(Occ0(c, BWT, i, true)!=cnt[4]){ //'T'
            cout<<"i:"<<i<<" c:"<<c<<" unEqual!"<<" Occ0:"<<Occ0(c, BWT, i, true)<<" Occ:"<<cnt[4]<<endl;
            getchar();
        }
    }
    cout<<"Occ pass!"<<endl;
}

void testBwtcpy(){
    uint l=131, r=137;
    uchar tbwt[r-l+1];
    // u32* ctbwt=new u32[(r-l+1)/16 +1];
    // memset(ctbwt,0,((r-l+1)/16 +1)*4);
    bwtcpy(tbwt, cBWT, l, r);
    // cout<<"tbwt:"<<ctbwt[0]<<endl;
    for(int i=0; i<r-l+1; i++){
        cout<<tbwt[i];
    }
    cout<<endl;
    for(uint i=l,j=0; i<=r; i++,j++){ //测试正确的
        cout<<BWT[i];
    }
    // delete[] ctbwt;
}

void test(){
    uchar c; bool* flag; bool* f[4]; uint cnt[5]; int nn[4];
    uint i=2741956596;
    // Occ(cBWT, i, cnt, 0, nn, flag, f, c);
    Occ(i, cnt, 0, nn, flag, f, c);
    i=C1[c+1] +cnt[c+1] -1;
    cout<<"i:"<<i<<endl;
    cout<<"cntA:"<<cnt[1]<<endl;
    cout<<"cntC:"<<cnt[2]<<endl;
    cout<<"cntG:"<<cnt[3]<<endl;
    cout<<"cntT:"<<cnt[4]<<endl;
}

int main0(int args, char** argv){
    coding(); //编码表 (这个一定要在最前面，如果没编码就会出现错误)
    // uint N1 =readHg19();
    // cout<<"finish readHg19(). N:"<<N<<endl; getchar();

    double time=0.0;
    srand48(9); //固定种子
    string file ="glocate.idx.hg19.32.256.cross";
    loadIdx(file); //直接载入索引
    // getchar(); //查看索引大小
    
    clock_t start, end; 
    int sz=0; uint* res0; uint* res1; uint* res2;
    // string fname="rpRead.txt"; //查询模式
    // string fname="rpRead.txt.2000-3000";
    // string fname="rpRead.txt.3000-4000";
    // string fname="rpRead.txt.4000-5000";
    string fname="rpRead.txt.2000-3000.kmer";
    // string fname="rpRead.txt.3000-4000.kmer";
    // string fname="rpRead.txt.4000-5000.kmer";
    nrd2=nrd;

    /* // genRef();                           //1.产生ref序列（随机
    N1 =readHg19();                           //2.给seq中赋值（真实
    cout<<"N1:"<<N1<<endl;
    genC_BWT_sSA();                           //产生C表、BWT串、SA采样数组、并存磁盘
    // int nSeeds =genRead_Seeds();           //产生seeds集合
    // test(); getchar();
    // testOcc();
    // testBwtcpy();
    cout<<"construct index finish!"<<endl;
    getchar(); */

    /* nSeeds =genRead_Seeds_Mmap(nrd);       //生成符合条件的reads
    L =new uint[nSeeds];
    R =new uint[nSeeds];
    // m =new uint[2*nSeeds];
    m =new uint[nSeeds];                      //只存累积量  
    // cout<<"1.nSeeds:"<<nSeeds<<endl;                 
    ofstream out(fname.c_str());              //将符合条件的reads写入文件
    for(uint i=0; i<nrd; i++){
        for(uint j=0; j<nSeeds; j++){
            out<<i<<"."<<j<<" "<<readSet[i]->seeds[j]<<" "<<readSet[i]->L[j]<<" "<<readSet[i]->R[j]<<" "<<readSet[i]->R[j]-readSet[i]->L[j]+1<<endl;
        }
        out<<endl;
    }
    out.close(); */

    // for(uint i=0; i<nrd; i++){             //count查询
    //     // if(readSet[i]->rep){
    //     //     cout<<"reads# "<<i<<" pos:"<<readSet[i]->pos<<endl;
    //     // }else{
    //     //     cout<<"reads. "<<i<<" pos:"<<readSet[i]->pos<<endl;
    //     // }
    //     for(int j=0; j<nSeeds; j++){
    //         count(readSet[i]->seeds[j], sLen, readSet[i]->L[j], readSet[i]->R[j], false);     
    //     }
    // }
    // cout<<endl;
    // getchar();

    ifstream in(fname.c_str());              //读入reads信息
    string data;
    int i=0, j=0, lines=0;
    // nSeeds=(rLen-sLen)/skip + ((((rLen-sLen)%skip)==0)?0:1); //有余数则加一，无余数则加零
    while(getline(in, data)){
        if(data==""){
            break;
        }
        lines++;
    }
    nSeeds=lines;                            //文件中空行前为seed个数
    // cout<<"2.nSeeds:"<<nSeeds<<endl;
    // getchar();
    for(int i=0; i<nrd; i++){
        readSet[i] =new struct read(nSeeds); //给readSet分配空间
    }
    in.clear();
    in.seekg(0,ios::beg);       //文件回到开头位置
    while(getline(in, data)){
        string *ss;
        if(data==""){           //一个read中含的几个seed结束
            i++;
            // nSeeds=j;        //根据读入文件，确定seed个数
            j=0;
            continue;
        }
        ss=covInfo(data,5);
        readSet[i]->L[j]=(uint)(atol(ss[2].c_str()));
        readSet[i]->R[j]=(uint)(atol(ss[3].c_str()));
        j++;
    }
    nrd2=i;                     //根据读入文件，确定read个数
    L =new uint[nSeeds];
    R =new uint[nSeeds];
    m =new uint[nSeeds];        //只存累积量  
    in.close();

    /* int diff=0;
    for(uint i=0; i<nrd2; i++){               //0.验证正确性
        diff=0;
        // L=readSet[i]->L;
        // R=readSet[i]->R;
        // memcpy(L,readSet[i]->L,sizeof(uint)*nrd2);
        // memcpy(R,readSet[i]->R,sizeof(uint)*nrd2);
        for(uint j=0; j<nSeeds; j++){
            L[j] =readSet[i]->L[j];
            R[j] =readSet[i]->R[j];
        }
        // res0 =locate(readSet[i]->L, readSet[i]->R, nSeeds, sz);    //locate查询
        // // res1 =glocate(L, R, nSeeds);       //glocate查询
        // res2 =glocatePlus(readSet[i]->L, readSet[i]->R, nSeeds, res0);   //glocate+查询
        res0 =locate(L, R, nSeeds, sz);    //locate查询
        res1 =glocate(L, R, nSeeds);       //glocate查询
        res2 =glocatePlus(L, R, nSeeds);   //glocate+查询
        // res2 =glocatePlus(L, R, nSeeds, res0);   //glocate+查询 (只能用于单一res0检查)
        for(uint j=0; j<sz; j++){
            if(res2[j]!=res0[j]){
                diff++;
                cout<<"第 "<<i<<" 个读长, 元素 "<<j<<" 结果不一致! ("<<res0[j]<<", "<<res2[j]<<")"<<endl;
                // getchar();
            }else{
                // cout<<j<<" 结果一致, "<<res2[j]<<endl;
            }
        }
        if(diff==0){
            cout<<i<<" Pass! 结果完全一致. (共"<<sz<<"个位置)"<<endl;
        }else{
            cout<<"    Error! 第 "<<i<<" 个读长, 有"<<diff<<"个位置不同, 共"<<sz<<"个位置"<<endl;
        }
    }
    delete[] res0;
    delete res1;
    delete[] res2; */
    // getchar();
    cout<<endl;
    
    start=clock();
    for(uint i=0; i<nrd2; i++){
        L=readSet[i]->L;
        R=readSet[i]->R;
        res0 =locate(L, R, nSeeds, sz);     //1.locate查询
    }
    end=clock();
    time=(double)(end-start)/CLOCKS_PER_SEC;
    cout<<"locate  nLF:"<<setw(8)<<nLF<<"    sz:"<<setw(7)<<sz<<"   time:"<<time<<endl;
    delete[] res0;
    // getchar();
    
    start=clock();
    for(uint i=0; i<nrd2; i++){
        L=readSet[i]->L;
        R=readSet[i]->R;
        res1 =glocate(L, R, nSeeds);        //2.glocate查询
    }
    end=clock();
    time=(double)(end-start)/CLOCKS_PER_SEC;
    cout<<"glocate nLF:"<<setw(8)<<nLF<<"   nLF1:"<<setw(6)<<nLF1<<"   time:"<<time<<endl;
    delete[] res1;
    // getchar();

    start=clock();
    for(uint i=0; i<nrd2; i++){
        L=readSet[i]->L;
        R=readSet[i]->R;
        res2 =glocatePlus(L, R, nSeeds);    //3.glocate+查询
    }
    end=clock();
    time=(double)(end-start)/CLOCKS_PER_SEC;
    cout<<"gloca++ nLF:"<<setw(8)<<nLF<<"   nLF1:"<<setw(6)<<nLF1<<"   time:"<<time<<endl;
    delete[] res2;

    return 1;
}

int main(int args, char** argv){
    coding(); //编码表 (这个一定要在最前面，如果没编码就会出现错误)
    // uint N1 =readHg19();
    // cout<<"finish readHg19(). N:"<<N<<endl; getchar();

    int task=0;
    cout<<"please input task (number:1~3):"<<endl;
    cout<<"1 index building"<<endl;
    cout<<"2 pattern generating"<<endl;
    cout<<"3 pattern locating"<<endl;
    cin>>task;
    if(task==1){
        cout<<"input genome file name:"<<endl;
        string file="";
        cin>>file;
        N1 =readHg19(file);                       //2.给seq中赋值（真实
        cout<<"N1:"<<N1<<endl;
        genC_BWT_sSA();                           //产生C表、BWT串、SA采样数组、并存磁盘
        cout<<"construct index finish!"<<endl;
    }else if(task==2){
        cout<<"input pattern set name:"<<endl;
        string fname="";
        cin>>fname;
        nSeeds =genRead_Seeds_Mmap(nrd);          //生成符合条件的reads
        L =new uint[nSeeds];
        R =new uint[nSeeds];
        // m =new uint[2*nSeeds];
        m =new uint[nSeeds];                      //只存累积量  
        // cout<<"1.nSeeds:"<<nSeeds<<endl;                 
        ofstream out(fname.c_str());              //将符合条件的reads写入文件
        for(uint i=0; i<nrd; i++){
            for(uint j=0; j<nSeeds; j++){
                out<<i<<"."<<j<<" "<<readSet[i]->seeds[j]<<" "<<readSet[i]->L[j]<<" "<<readSet[i]->R[j]<<" "<<readSet[i]->R[j]-readSet[i]->L[j]+1<<endl;
            }
            out<<endl;
        }
        out.close();
        cout<<"generate pattern finish!"<<endl;
    }else if(task==3){
        cout<<"input pattern set name:"<<endl;
        string fname="";
        cin>>fname;
        cout<<"input index file name:"<<endl;
        string findex="";
        cin>>findex;
        double time=0.0;
        srand48(9); //固定种子
        loadIdx(findex); //直接载入索引
        clock_t start, end; 
        int sz=0; uint* res0; uint* res1; uint* res2;
        // string fname="rpRead.txt"; //查询模式
        // string fname="rpRead.txt.2000-3000";
        // string fname="rpRead.txt.3000-4000";
        // string fname="rpRead.txt.4000-5000";
        // string fname="rpRead.txt.2000-3000.kmer";
        // string fname="rpRead.txt.3000-4000.kmer";
        // string fname="rpRead.txt.4000-5000.kmer";
        nrd2=nrd;

        ifstream in(fname.c_str());              //读入reads信息
        string data;
        int i=0, j=0, lines=0;
        // nSeeds=(rLen-sLen)/skip + ((((rLen-sLen)%skip)==0)?0:1); //有余数则加一，无余数则加零
        while(getline(in, data)){
            if(data==""){
                break;
            }
            lines++;
        }
        nSeeds=lines;                            //文件中空行前为seed个数
        // cout<<"2.nSeeds:"<<nSeeds<<endl;
        // getchar();
        for(int i=0; i<nrd; i++){
            readSet[i] =new struct read(nSeeds); //给readSet分配空间
        }
        in.clear();
        in.seekg(0,ios::beg);       //文件回到开头位置
        while(getline(in, data)){
            string *ss;
            if(data==""){           //一个read中含的几个seed结束
                i++;
                // nSeeds=j;        //根据读入文件，确定seed个数
                j=0;
                continue;
            }
            ss=covInfo(data,5);
            readSet[i]->L[j]=(uint)(atol(ss[2].c_str()));
            readSet[i]->R[j]=(uint)(atol(ss[3].c_str()));
            j++;
        }
        nrd2=i;                     //根据读入文件，确定read个数
        L =new uint[nSeeds];
        R =new uint[nSeeds];
        m =new uint[nSeeds];        //只存累积量  
        in.close();

        /* int diff=0;
        for(uint i=0; i<nrd2; i++){               //0.验证正确性
            diff=0;
            // L=readSet[i]->L;
            // R=readSet[i]->R;
            // memcpy(L,readSet[i]->L,sizeof(uint)*nrd2);
            // memcpy(R,readSet[i]->R,sizeof(uint)*nrd2);
            for(uint j=0; j<nSeeds; j++){
                L[j] =readSet[i]->L[j];
                R[j] =readSet[i]->R[j];
            }
            // res0 =locate(readSet[i]->L, readSet[i]->R, nSeeds, sz);    //locate查询
            // // res1 =glocate(L, R, nSeeds);       //glocate查询
            // res2 =glocatePlus(readSet[i]->L, readSet[i]->R, nSeeds, res0);   //glocate+查询
            res0 =locate(L, R, nSeeds, sz);    //locate查询
            res1 =glocate(L, R, nSeeds);       //glocate查询
            res2 =glocatePlus(L, R, nSeeds);   //glocate+查询
            // res2 =glocatePlus(L, R, nSeeds, res0);   //glocate+查询 (只能用于单一res0检查)
            for(uint j=0; j<sz; j++){
                if(res2[j]!=res0[j]){
                    diff++;
                    cout<<"第 "<<i<<" 个读长, 元素 "<<j<<" 结果不一致! ("<<res0[j]<<", "<<res2[j]<<")"<<endl;
                    // getchar();
                }else{
                    // cout<<j<<" 结果一致, "<<res2[j]<<endl;
                }
            }
            if(diff==0){
                cout<<i<<" Pass! 结果完全一致. (共"<<sz<<"个位置)"<<endl;
            }else{
                cout<<"    Error! 第 "<<i<<" 个读长, 有"<<diff<<"个位置不同, 共"<<sz<<"个位置"<<endl;
            }
        }
        delete[] res0;
        delete res1;
        delete[] res2; */
        // getchar();
        cout<<endl;
        
        start=clock();
        for(uint i=0; i<nrd2; i++){
            L=readSet[i]->L;
            R=readSet[i]->R;
            res0 =locate(L, R, nSeeds, sz);     //1.locate查询
        }
        end=clock();
        time=(double)(end-start)/CLOCKS_PER_SEC;
        cout<<"locate  nLF:"<<setw(8)<<nLF<<"    sz:"<<setw(7)<<sz<<"   time:"<<time<<endl;
        delete[] res0;
        // getchar();
        
        start=clock();
        for(uint i=0; i<nrd2; i++){
            L=readSet[i]->L;
            R=readSet[i]->R;
            res1 =glocate(L, R, nSeeds);        //2.glocate查询
        }
        end=clock();
        time=(double)(end-start)/CLOCKS_PER_SEC;
        cout<<"glocate nLF:"<<setw(8)<<nLF<<"   nLF1:"<<setw(6)<<nLF1<<"   time:"<<time<<endl;
        delete[] res1;
        // getchar();

        start=clock();
        for(uint i=0; i<nrd2; i++){
            L=readSet[i]->L;
            R=readSet[i]->R;
            res2 =glocatePlus(L, R, nSeeds);    //3.glocate+查询
        }
        end=clock();
        time=(double)(end-start)/CLOCKS_PER_SEC;
        cout<<"gloca++ nLF:"<<setw(8)<<nLF<<"   nLF1:"<<setw(6)<<nLF1<<"   time:"<<time<<endl;
        delete[] res2;
    }else{
        cout<<"Unavailable options, Effloc will exit!"<<endl;
    }

    return 1;
}