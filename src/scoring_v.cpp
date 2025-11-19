/////////////////////

//      lastest scoring function        
/*  use hashing method to speed up input scoring */
/*  2 mismatch                                   */
/*  new background scoring                       */
/*  allow multiple bases for gaps                */

//////////////////////


#include<iostream>
#include<math.h>
#include<algorithm>
#include <stdlib.h>
#include<vector>
#include<cstring>
#include<omp.h>

using namespace std;

#define max_numseq 1001
#define max_seqlen 5001
#define max_nummof 200000

char filefreq6[100]="AG.6.freq";
char filefreq8[100]="AG.8.freq";
char fileinput[100]="d1.fasta";
char filemotif[100]="d1.dat";

#define size 5
#define max 1024  /*4^5*/
#define nummism 2

int num;
char seq[max_numseq][max_seqlen];
char gmotif[100];

double table6[5000];
double tab6[20000];
double table8[70000];
double tab8[400000];
double eps;

int h6=12,h8=16;

int mark[max_numseq][max_seqlen];
int sublen;

char fileconfig[100];

int legal(char c){
    return ('0'<=c&&c<='9')||('a'<=c&&c<='z')||('A'<=c&&c<='Z');
}

void trim(char *s){
    int m=strlen(s);
    while(m&&!legal(s[m-1])) s[--m]=0;
}

void readParameter(){
    FILE *f=fopen(fileconfig,"r");
    char s[100];
    fgets(s,100,f);sscanf(s,"%i",&sublen);    
    fgets(s,100,f);fgets(s,100,f);fgets(s,100,f);
    fgets(filefreq6,100,f);
    trim(filefreq6);
    fgets(filefreq8,100,f);
    trim(filefreq8);    
    fgets(fileinput,100,f);
    trim(fileinput);    
    fgets(filemotif,100,f);
    trim(filemotif);    
    fclose(f);
}

bool junk(char c){
    return !('a'<=c&&c<='z')&&!('A'<=c&&c<='Z');
}

void capitalize(char *s){
    int i,t=strlen(s);
    for (i=0;i<t;i++)
        if ('a'<=s[i]&&s[i]<='z') s[i]=s[i]-'a'+'A';
}

void readData(){
    FILE *f;
    num=0;
       
/*    while(gets(s)!=NULL){
        gets(seq[num]);
        cout<<strlen(seq[num])<<endl;
        num++;
    }    */
    f=fopen(fileinput,"r");
    num=0;
    char s[200];
    while(fgets(s,max_seqlen,f)!=NULL&&fgets(seq[num],max_seqlen,f)!=NULL){
        if (junk(seq[num][strlen(seq[num])-1])) seq[num][strlen(seq[num])-1]=0;
        capitalize(seq[num]);        
//        cout<<strlen(seq[num])<<endl;
        num++;
    }
}

int vt(char c){
    switch (c){
        case 'A':return 0;
        case 'C':return 1;
        case 'G':return 2;
    }
    return 3;
}

int con(char *s){
    int i,t=0,k=1;
    for (i=0;i<size;i++){
        t+=vt(s[size-i-1])*k;k*=4;
    }
    return t;
}

typedef struct{
    int se;int po;int mism;
} posi;

vector<posi> a[max];
int b[max];
char s[size+1];
char base[]="ACGT";

// # of mismatched between s1 and s2, not more than nummism
int check(char *s1,char *s2){
    int i,t=0;
    int m=strlen(s1);
    for (i=0;i<m;i++)
        if (s1[i]!=s2[i]){
            t++;
            if (t>nummism) return -1;
        }
    return t;
}

void tryo(int x){
    int i,j;
    if (x==size){
        s[x]=0;
        int y=con(s);
        for (i=0;i<num;i++){
            const int m=strlen(seq[i]);  // Cache strlen
            for (j=0;j<m-size+1;j++){
                int t=check(s,seq[i]+j);
                if (t!=-1){
                    posi pp;
                    pp.se=i;pp.po=j;pp.mism=t;
//                    a[y][b[y]].se=i;a[y][b[y]].po=j;a[y][b[y]].mism=t;
                    a[y].push_back(pp);
                    b[y]++;
                }
            }
        }
        return;
    }
    for (i=0;i<4;i++){
        s[x]=base[i];tryo(x+1);
    }
}

int stnum;
#pragma omp threadprivate(mark, stnum)

void process(){
    int i,j;
    for (i=0;i<max;i++){
        b[i]=0;
        a[i].reserve(200);  // Reserve capacity to avoid reallocations
    }

    tryo(0);  // construct the hashtable

    for (i=0;i<num;i++){
        const int m=strlen(seq[i]);  // Cache strlen
        for (j=0;j<m;j++) mark[i][j]=0;
    }
    stnum=0;
}

//read frequency files
void readFreqFile(){
    int i;
    double to6=0,to8=0;
    
    FILE *f=fopen(filefreq6,"r");
    for (i=0;i<(1<<12);i++){
        char s1[100];
        int x;
        fscanf(f,"%s %f",s1,&x);
        table6[i]=x;        
        to6+=table6[i];
    }
    fclose(f);

    f=fopen(filefreq8,"r");
    for (i=0;i<(1<<16);i++){
        char s1[100];
        int x;
        fscanf(f,"%s %i\n",s1,&x);
        table8[i]=x;
        to8+=table8[i];
    }   
    fclose(f);
    
    for (i=0;i<(1<<12);i++) table6[i]/=to6;
    for (i=0;i<(1<<16);i++) table8[i]/=to8;
    
    eps=1.0/to8;
}

//calculate the number of matched between s1 and prefix of s2
int count1(char *s1,char *s2,int th){
    int i,l=0;
    int n=strlen(s1);
    for (i=0;i<n;i++)
        if (s1[i]!='N'&&s1[i]!=s2[i]){
            if (l==th) return -1;            
            l++;
        }
    return l;
}

int best,nbest=1;

int base_hash(char c){
    switch(c){
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: return 4;
    }
}

//transform DNA sequence into 5th ary
int trans5(char *s,int h){
    int i,mu4=1;
    int t=0;
    for (i=h-1;i>=0;i--){
        t+=mu4*base_hash(s[i]);
        mu4*=5;
    }
    return t;
}

//transform DNA sequence into 4th ary
int trans(char *s,int h){
    int i,mu4=1;
    int t=0;
    for (i=h-1;i>=0;i--){
        t+=mu4*base_hash(s[i]);
        mu4*=4;
    }
    return t;
}

double cal(char *s){
    int i,h=6;

    int m=strlen(s);
    
    if (strlen(s)>=8) h=8;    
    int win=trans5(s,h);
    double exp=tab8[win];
    if (h==6) exp=tab6[win];

    int mu5=1;
    for (i=0;i<h-1;i++) mu5*=5;
    
    for (i=h;i<m;i++){
        win=(win-base_hash(s[i-h]))/5;win=win+base_hash(s[i])*mu5;

        double t1=tab8[win];
        if (h==6) t1=tab6[win];

        double t2=0;

            char tm=s[i];
            s[i]='N';
            if (h==8) t2+=tab8[win+(4-base_hash(s[i]))*mu5];
            if (h==6) t2+=tab6[win+(4-base_hash(s[i]))*mu5];
            s[i]=tm;

        if (t2>0) exp*=t1/t2;
    }        
    return exp;
}

void permut(int x,char *s){
    char base[]="ACGTN";
    int i,j;
    if (x==6){
        int x=trans5(s,6);
        int ok=0;
        for (i=5;i>=0;i--)
            if (s[i]=='N'){
                ok=1;break;
            }
        if (ok){
            tab6[x]=0;
            for (j=0;j<4;j++){
                s[i]=base[j];
                tab6[x]+=tab6[trans5(s,6)];          
            }
        }else tab6[x]=table6[trans(s,6)];
    }

    if (x==8){
        int x=trans5(s,8);
        int ok=0;
        for (i=7;i>=0;i--)
            if (s[i]=='N'){
                ok=1;break;
            }
        if (ok){
            tab8[x]=0;
            for (j=0;j<4;j++){
                s[i]=base[j];
                tab8[x]+=tab8[trans5(s,8)];
            }
        }else tab8[x]=table8[trans(s,8)];
        return;
    }
    
    if (x<8)
        for (i=0;i<5;i++){
            s[x]=base[i];
            permut(x+1,s);
        }
}

void preprocess(){
    char s[100];
    permut(0,s);
}

//calculate probability of gapped motif appears in background sequence
//beta

int counter=0;

typedef struct{
    int min;int dd;
} str;

str pers[max_numseq];

// find all occurences of s bases on hashing value at position
void find(int pos,char *s){
    int i,j;
    int fir=con(s+pos);
    int m=strlen(s);
    
    for (i=0;i<b[fir];i++){
        int x=a[fir][i].se;
        int y=a[fir][i].po;
        int numm=a[fir][i].mism;
        
        if (y-pos<0) continue;        
        if (mark[x][y-pos]==stnum) continue;
        
        if (pers[x].min<numm) continue;
        
        mark[x][y-pos]=stnum;
        
        for (j=y-pos;j<y&&numm<=pers[x].min;j++)
            if (s[j-y+pos]!='N'&&seq[x][j]!=s[j-y+pos]) numm++;
               
        for (j=y+size;j<y+m-pos&&numm<=pers[x].min;j++)
            if (s[j-y+pos]!='N'&&seq[x][j]!=s[j-y+pos]) numm++;
        
        if (numm==pers[x].min) pers[x].dd++;
        if (numm<pers[x].min){
            pers[x].dd=1;pers[x].min=numm;
        }
    }    
}

// brute force
void brutef(int x,char *s,int min){
    int L=strlen(seq[x]);
    int m=strlen(s);

    int cc=0;
    
    for (int j=0;j<L-m+1;j++){
        int k=count1(s,seq[x]+j,min);
        if (k==min) cc++;
        if (k>-1&&k<min){
            min=k;cc=1;
        }
    }    
    pers[x].min=min;pers[x].dd=cc;
}

// probability of s in bg with e mismatches
double background_pro(char *s,int e,int x){
    int m=strlen(s);

    if (x==m) return cal(s);
    
    if (s[x]=='N') return background_pro(s,e,x+1);
    else{
        double to=0;
        if (e){
            char tg=s[x];
            s[x]='N';
            to+=background_pro(s,e-1,x+1);                                
            s[x]=tg;
        }
        return to+background_pro(s,e,x+1);
    }
}

double sigma(char *s){
    int i;
    double re=0;    
    for (i=0;i<num;i++){
        double tmp=background_pro(s,pers[i].min,0);
        if (tmp<eps) tmp=eps;
        re+=log(1.0/(tmp*strlen(seq[i])));
    }
    return re;
}

double beta(char *s){
    int i;
    int m=strlen(s);
    
    best=m;nbest=1;
    int beg[2];
    int st=0,dd=0;
    
    for (i=0;i<m;i++)
        if (s[i]!='N'){
            dd++;
            if (dd==size){
                beg[st]=i-size+1;dd=0;st++;
                if (st==2) break;
            }                 
        }else dd=0;
        
    stnum++;
    for (i=0;i<num;i++) pers[i].min=1000000;

    if (st==2){
        find(beg[0],s);find(beg[1],s);
    }
    
    for (i=0;i<num;i++)
        if (pers[i].min>5) brutef(i,s,pers[i].min);

    best=m;nbest=1;    

    for (i=0;i<num;i++){
        if (best==pers[i].min) nbest+=pers[i].dd;
        if (best>pers[i].min){
            best=pers[i].min;nbest=pers[i].dd;
        }
    }

    double sumlen=0;
    // Cache strlen to avoid repeated calls
    for (i=0;i<num;i++) sumlen+=strlen(seq[i]);
    
        double tmp=background_pro(s,best,0);
        if (tmp<eps) tmp=eps;
    
    return log(nbest/(tmp*sumlen));
}

char motif[max_nummof][30];
int nummof;
double score[max_nummof];
int tt[max_nummof];

void qsort(int l,int r){
    int i=l,j=r;
    double w=score[(i+j)/2];
    do{
        while(score[i]>w) i++;
        while(score[j]<w) j--;
        if (i<=j){
            int tg=tt[i];tt[i]=tt[j];tt[j]=tg;
            double tg2=score[i];score[i]=score[j];score[j]=tg2;
            i++;j--;
        }
    }while(i<=j);
    if (i<r) qsort(i,r);
    if (l<j) qsort(l,j);
}

typedef struct{
    int seq;int pos;
} subpos;
subpos c[max_nummof][max_numseq];

    int map(char c){
        switch (c){
            case 'A':return 1;
            case 'C':return 2;     
            case 'G':return 4;
        }
        return 8;
    }

int numsub[max_nummof];

void cluster(int x){
    int i,j;
    int m=strlen(motif[x]);
    for (i=0;i<m;i++)
        if (motif[x][i]=='N'){
            int dd[300];
            dd[(int)'A']=0;dd[(int)'C']=0;dd[(int)'G']=0;dd[(int)'T']=0;
            for (j=0;j<numsub[x];j++) dd[(int)seq[c[x][j].seq][c[x][j].pos+i]]=1;
            motif[x][i]=dd[(int)'A']*map('A')+dd[(int)'C']*map('C')+dd[(int)'G']*map('G')+dd[(int)'T']*map('T');
        }
}

void print(char *s){
    int m=strlen(s),i;
    for (i=0;i<m;i++)
        if (s[i]>16) cout<<s[i];
        else{
            cout<<'[';
            int x=(int)s[i];
            if ((1&(x>>0))==1) cout<<'A';
            if ((1&(x>>1))==1) cout<<'C';
            if ((1&(x>>2))==1) cout<<'G';
            if ((1&(x>>3))==1) cout<<'T';
            cout<<']';
        }
}

void ssort(int l,int r){
    int i=l,j=r;
    char w[max_seqlen];
    strcpy(w,motif[tt[(i+j)/2]]);
    
    do{
        while(strcmp(motif[tt[i]],w)<0) i++;
        while(strcmp(w,motif[tt[j]])<0) j--;        
        if (i<=j){
            int tg=tt[i];tt[i]=tt[j];tt[j]=tg;
            i++;j--;
        }
    }while(i<j);
    if (i<r) ssort(i,r);
    if (l<j) ssort(l,j);
}   
    
int main(int narg,char *arg[]){
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    time_t seconds,end;
    seconds = time (NULL);

    if (narg<2){
        cout<<"you need to specify config file";return 0;
    }
    strcpy(fileconfig,arg[1]);
    
    readParameter();
    readData();
    readFreqFile();

    FILE *f=fopen(filemotif,"r");
    nummof=0;                
    
    while(fscanf(f,"%s",motif[nummof])!=EOF){
        fscanf(f,"%i",&numsub[nummof]);
        if (strcmp(motif[nummof],"end")==0) break;
        for (int i=0;i<numsub[nummof];i++) fscanf(f,"%i %i",&c[nummof][i].seq,&c[nummof][i].pos);
        if (nummof==0||strcmp(motif[nummof],motif[nummof-1])!=0) nummof++;
    }

    int i;
    for (i=0;i<nummof;i++) tt[i]=i;        
    ssort(0,nummof-1);
    
    bool del[max_nummof];
    memset(del,false,nummof);
    for (i=1;i<nummof;i++)
        if (strcmp(motif[tt[i]],motif[tt[i-1]])==0) del[tt[i]]=true;
        
    process();
    preprocess();

    //scoring - parallelized for multi-core performance
    #pragma omp parallel for schedule(dynamic) if(nummof>10)
    for (i=0;i<nummof;i++)
        if (!del[i]){
//            cout<<motif[i]<<endl;
            score[i]=beta(motif[i])+sigma(motif[i]);
        }else score[i]=-10000000;        

    for (i=0;i<nummof;i++) tt[i]=i;
    
    qsort(0,nummof-1);
    while(nummof&&score[nummof-1]==-10000000) nummof--;
    
    //clustering
    for (i=0;i<nummof;i++) cluster(tt[i]);
    
    for (i=0;i<nummof;i++){
        int x=tt[i];
        cout<<"MOTIF : ";
        print(motif[x]);cout<<" "<<score[i]<<endl;
        int len=strlen(motif[x]);
        for (int j=0;j<numsub[x];j++){
            int y=c[x][j].seq;
//            cout<<strlen(seq[y])<<endl;
            cout<<y<<",-"<<strlen(seq[y])-c[x][j].pos-1<<",";
            int l;
            for (l=0;l<len;l++) cout<<seq[y][c[x][j].pos+l];
            cout<<endl;
        }
            cout<<"========"<<endl;
    }
    
    end = time (NULL);
    cout<<"Total time: "<<end-seconds<<endl;
}
