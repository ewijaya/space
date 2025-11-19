/////////////////////

//  lastest program
//  extend the size of hash tables
//  use vectors to make program flexible

/////////////////////

#include<iostream>
#include<math.h>
#include<vector>
#include<cstring>

using namespace std;

#define max_numseq 1000
#define max_seqlen 5001
#define max_sublen 5
#define num_string 4400
#define d 1

int sublen=4;
int coverage=10;
double q=0.5;
int motiflen=15;

char base[4]={'A','G','C','T'};
char seq[max_numseq][max_seqlen];
char entry[max_sublen+1];
int num_entry;
typedef struct{
    int seq;int pos;
} Ipair;

vector<Ipair> list1[num_string];
vector<Ipair> list2[num_string];
vector<Ipair> list3[num_string];

int dd1[num_string],dd2[num_string],dd3[num_string];
int num;
int lp=15;
int pos1[max_numseq],pos2[max_numseq];
int qu;

///////////////////////////////////////////////////////////////

vector<int> a[max_numseq][max_seqlen];    // seq number
vector<int> b[max_numseq][max_seqlen];    // possition
int c[max_numseq][max_seqlen];

char fileconfig[100];

void read_parameter(){
    FILE *f=fopen(fileconfig,"r");
    
    fscanf(f,"%i %i %lf %i",&sublen,&coverage,&q,&motiflen);
    fclose(f);
}

int check(int can,int x,int y,int k){
    int i;
    for (i=0;i<c[can][x];i++)
        if (a[can][x][i]==y&&b[can][x][i]==k) return 0;
        
/*    int jj=a[can][x][mx3];
    while(jj!=-1){
        for (i=0;i<c2[jj];i++)
            if (a2[jj][i]==y&&b2[jj][i]==k) return 0;        
        jj=a2[jj][mx3];
    }*/
    return 1;
}

int total=0;
int jump=0;

void push_back1(int can,int x,int y,int k){
    if (check(can,x,y,k)){
/*        if (c[can][x]==mx3){
              int jj=a[can][x][mx3];
              if (jj==-1){
                    assert(jump<mx4);

                    a[can][x][mx3]=jump;jj=jump;jump++;
                    if (jump<mx4) a2[jump][mx3]=-1;
                    c2[jj]=1;a2[jj][0]=y;b2[jj][0]=k;                    
              }else{
                    while(a2[jj][mx3]!=-1) jj=a2[jj][mx3];
                    if (c2[jj]<mx3){
                        a2[jj][c2[jj]]=y;b2[jj][c2[jj]]=k;
                        c2[jj]++;
                    }else{           
                        assert(jump<mx4);             

                        a2[jj][mx3]=jump;
                        c2[jump]=1;a2[jump][0]=y;b2[jump][0]=k;
                                                
                        jump++;
                        if (jump<mx4) a2[jump][mx3]=-1;                        
                    }                                                                   
              }
        }else*/{
            a[can][x].push_back(y);
            b[can][x].push_back(k);
            c[can][x]++;
        }
        total++;
    }    
}

void push_back(int can,int x,int y,int k){
    push_back1(can,x,y,k);push_back1(y,k,can,x);
}

void pre_align(){
    int can,i;
    for (can=0;can<num;can++)
        for (i=0;i<max_seqlen;i++){
            c[can][i]=0;//a[can][i][mx3]=-1;
        }
    jump=0;
}

void alignment(int x,int y,int z,int x1,int y1,int z1){
    int k;               
                int t1=strlen(seq[x]);
                int t2=strlen(seq[x1]);
                                
                for (k=0;k<=motiflen-(z-y+sublen);k++)
                    if (y>=k&&y1>=k)
                    if (y-k+motiflen<=t1&&y1-k+motiflen<=t2)
                    push_back(x,y-k,x1,y1-k);
}

//////////////////////////////////////////////////////////////////////

int compare(char *s1,char *s2,int k){
    int i,dd=0;
    for (i=0;i<k;i++){
        if (s1[i]!=s2[i]) dd++;
        if (dd>d) return 0;
    }
    return 1;
}

void find(int x){
    int i,j;
    dd1[x]=0;
        
    for (i=0;i<num;i++){
        int m=strlen(seq[i]);
        for (j=0;j<m;j++)
            if (strncmp(entry,seq[i]+j,sublen)==0){
                Ipair ii;
                ii.seq=i;ii.pos=j;
                list1[x].push_back(ii);
                dd1[x]++;
            }
    }

    for (i=0;i<num;i++){
        int m=strlen(seq[i]);
        for (j=0;j<m-sublen+1;j++)
            if (compare(entry,seq[i]+j,sublen)){
                Ipair ii;
                ii.seq=i;ii.pos=j;
                list2[x].push_back(ii);
                
//                list2[x][dd2[x]][0]=i;list2[x][dd2[x]][1]=j;
                dd2[x]++;
            }
    }
}

int vt(char c){
    switch (c){
        case 'A':return 0;
        case 'G':return 1;
        case 'C':return 2;
    }
    return 3;
}

int con(char *s){
    int i,t=0,k=1;
    for (i=0;i<sublen;i++){
        t+=vt(s[sublen-i-1])*k;k*=4;
    }
    return t;
}

void find_hashtable1(){
    int i,j,l,g;

    char s[sublen+1];
    
    for (i=0;i<num;i++)
        for (j=0;j<(int)strlen(seq[i])-sublen+1;j++){
            strncpy(s,seq[i]+j,sublen);
            s[sublen]=0;
            
            int x=con(s);
                Ipair ii;
                ii.seq=i;ii.pos=j;
                list1[x].push_back(ii);
            
//            list1[x][dd1[x]][0]=i;list1[x][dd1[x]][1]=j;
                dd1[x]++;            

                list2[x].push_back(ii);
                list3[x].push_back(ii);
                    
//            list2[x][dd2[x]][0]=i;list2[x][dd2[x]][1]=j;
            dd2[x]++;
//            list3[x][dd3[x]][0]=i;list3[x][dd3[x]][1]=j;
            dd3[x]++;
            
            for (l=0;l<sublen;l++){
                char ch=s[l];
                for (g=0;g<4;g++)
                    if (base[g]!=ch){
                        s[l]=base[g];
                        x=con(s);

                        Ipair ii;
                        ii.seq=i;ii.pos=j;
                        list2[x].push_back(ii);                        
//                        list2[x][dd2[x]][0]=i;list2[x][dd2[x]][1]=j;
                        dd2[x]++;
                        
                        if (base[g]>ch){
                            list3[x].push_back(ii);
//                            list3[x][dd3[x]][0]=i;list3[x][dd3[x]][1]=j;
                            dd3[x]++; 
                        }
                    }
                s[l]=ch;
            }
        }
}

typedef struct{
    int seq;int st1;int st2;
} triple;

char entries[num_string][max_sublen+1];   // AAAAA to TTTTT

void find2(){
    int i,j,l,i1,j1;
    int x,y,z;
                                  
    int seqq,st1,st2;
    for (i=0;i<num_entry;i++)
        for (j=0;j<dd1[i];j++){
            x=list1[i][j].seq;
            y=list1[i][j].pos;

            for (l=sublen;l<=lp-sublen;l++)
                if (y+l+sublen-1<(int)strlen(seq[x])){
                    z=con(seq[x]+y+l);
                    
                    int ok=0;
                    j1=0;
                    for (i1=0;i1<dd3[i];i1++){
                        while(j1<dd2[z]&&list3[i][i1].seq>list2[z][j1].seq||
                        (list3[i][i1].seq==list2[z][j1].seq&&list3[i][i1].pos+l>list2[z][j1].pos)) j1++;
                        if (j1>=dd2[z]) break;

                        if (list3[i][i1].seq!=x||list3[i][i1].pos!=y)
                        if (list3[i][i1].seq==list2[z][j1].seq&&list3[i][i1].pos+l==list2[z][j1].pos){
                                if (ok==0){
                                    seqq=x;
                                    st1=y;
                                    st2=y+l;
                                    ok=1;                                    
                                }

                                alignment(seqq,st1,st2,list3[i][i1].seq,list3[i][i1].pos,list2[z][j1].pos);
                            }
                    }
                }
        }
}

void sort(int x,int y){
    int i,j,tg;
    for (i=0;i<c[x][y];i++)
        for (j=i+1;j<c[x][y];j++)
            if (a[x][y][i]>a[x][y][j]||(a[x][y][i]==a[x][y][j]&&b[x][y][i]>b[x][y][j])){
                tg=a[x][y][i];
                a[x][y][i]=a[x][y][j];
                a[x][y][j]=tg;
                tg=b[x][y][i];
                b[x][y][i]=b[x][y][j];
                b[x][y][j]=tg;
            }
}

#define mx3 500
int a3[mx3*10],b3[mx3*10];
int cnt;

void extract(int can,int x){
    int i;
    cnt=0;
    for (i=0;i<c[can][x];i++){
        a3[cnt]=a[can][x][i];
        b3[cnt]=b[can][x][i];
        cnt++;
    }
    
/*    int jj=a[can][x][mx3];    
    while(jj!=-1){
        for (i=0;i<c2[jj];i++){
            a3[cnt]=a2[jj][i];
            b3[cnt]=b2[jj][i];
            cnt++;
        }
        jj=a2[jj][mx3];
    }*/
}

void result2(){
    int i,j,l;
    int can=0;
    cout<<qu<<endl;
    int count=0;

    int longest=0;

    for (can=0;can<num;can++)
        for (i=0;i<=(int)strlen(seq[can])-motiflen;i++)
            if (c[can][i]){
                extract(can,i);                
                                
/*                int mk[max_numseq];
                for (j=0;j<num;j++) mk[j]=0;
                int kk=0;               
                for (j=0;j<cnt;j++){               
                    int x=a3[j];
                    if (mk[x]==0){
                        mk[x]=1;kk++;
                    }
                }
                if (kk<qu) continue;*/
                
                count++;
                for (j=i;j<i+motiflen;j++) cout<<seq[can][j];

                cout<<endl<<sublen<<" "<<coverage<<" "<<can<<" "<<i<<endl;

                int dd=0;
                int *e;
                e=(int*)malloc((motiflen+2)*sizeof(int));

                int *g;
                g=(int*)malloc((motiflen+2)*sizeof(int));

                int ttt=0;
                for (j=0;j<cnt;j++){
                    int x=a3[j];
                    int y=b3[j];
                    for (l=0;l<motiflen;l++)
                        if (seq[can][i+l]==seq[x][y+l]) e[l]=0;
                        else e[l]=1;

                    e[motiflen]=0;
                    for (l=motiflen-1;l>=0;l--) e[l]=e[l]+e[l+1];

                    int r=0;
                    dd=0;
                    for (l=0;l<=motiflen-sublen;l++)
                        if (e[l]-e[l+sublen]<2){
                            if (dd&&g[dd-1]+sublen-1>=l) r+=sublen-(g[dd-1]+sublen-l);
                            else r+=sublen;
                            g[dd]=l;dd++;
                        }
                    if (r>=coverage){
                        cout<<x<<" "<<y<<" ";
                      cout<<dd;
                      for (l=0;l<dd;l++) cout<<" "<<g[l];
                      cout<<endl;
                        ttt++;
                    }
                }
                free(e);
                free(g);
              cout<<0<<" "<<0<<" "<<-1<<endl<<endl<<endl;
                if (longest<ttt) longest=ttt;
            }            
    cout<<"end"<<endl;
    cout<<"# of motifs: "<<count<<endl;    
    cout<<"longest set: "<<longest<<endl;
}

void result(){
    int i,j;
    int can=0;
    
    for (can=0;can<num;can++)
        for (i=0;i<=(int)strlen(seq[can])-motiflen;i++)
            if (c[can][i]){
                cout<<"Motif candidate "<<can<<" "<<i<<" ";
                for (j=i;j<i+motiflen;j++) cout<<seq[can][j];
                cout<<endl;
                cout<<"No of instances: "<<c[can][i]<<endl;            
                sort(can,i);
                for (j=0;j<c[can][i];j++) cout<<a[can][i][j]<<" "<<b[can][i][j]<<endl;
                cout<<endl<<endl;
            }
}

void capitalize(char *s){
    int i,t=strlen(s);
    for (i=0;i<t;i++)
        if ('a'<=s[i]&&s[i]<='z') s[i]=s[i]-'a'+'A';
}

void readdata(){
    num=0;
    char s[200];
    while(gets(s)!=NULL){
        gets(seq[num]);
        capitalize(seq[num]);
//        cout<<strlen(seq[num])<<endl;                
        num++;
    }    
    qu=(int)ceil(num*q);
    qu=qu-1;
}

void sub_gen(int x){
    int i;
    if (x==sublen){
        strcpy(entries[num_entry],entry);num_entry++;
        return;
    }
    for (i=0;i<4;i++){
        entry[x]=base[i];
        sub_gen(x+1);
    }
}

int main(int narg,char *arg[]){
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    if (narg<2){
        cout<<"you need to specify configuration file";return 0;
    }
    strcpy(fileconfig,arg[1]);

    time_t seconds,end;
    seconds = time (NULL);
          
    read_parameter();
    readdata();
    
    num_entry=0;

    sub_gen(0);

    find_hashtable1();
    
    pre_align();    
    find2();

    result2();
    
    end = time (NULL);
    cout<<"Total time: "<<end-seconds<<endl;
}
