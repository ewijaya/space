/*
    Frequent pattern mining
    length of arrays is arbitrary
*/

#include<iostream>
#include<vector>
#include<cstring>
#include<cassert>

#define max 1000  //max length of item set
#define max2 1000000
#define max3 30  //max length of motif

using namespace std;

vector<int> set1[max2],set2[max2];
long dset1[max2],dset2[max2];
long dd;
int moflen,cov;

typedef struct{
    int seq;int pos;
} subpos;

//input is a set of sorted set of integers
long intersect(vector<int>& a,vector<int>& b,long m,long n,vector<int>& c){
    long i=0,j=0,dd=0;
    while(i<m&&j<n){
        if (a[i]==b[j]){
            c.push_back(a[i]);dd++;
            i++;j++;
        }else if (a[i]>b[j]) j++;else i++;
    }
    return dd;
}

long del[max2];
void merge(long a[][max3],long b[],long n,long q){
    long i,j;
    vector<int> c[max];
    long d[max];
    
    long mx=0,mn=max;
    
    for (i=0;i<n;i++){
        if (a[i][b[i]-1]>mx) mx=a[i][b[i]-1];
        if (a[i][0]<mn) mn=a[i][0];
    }
    
    for (j=0;j<=mx;j++) d[j]=0;
    
    for (i=0;i<n;i++)
        for (j=0;j<b[i];j++){
            long x=a[i][j];
            c[x].push_back(i);d[x]++;
        }

    for (i=0;i<dd;i++) set1[i].clear();
    for (i=0;i<dd;i++) set2[i].clear();
    
    dd=0;
    for (i=mx;i>=0;i--){
        if (d[i]<q) continue;
        
        long ok=1,ok2=1;
        
        int t=dd;
        for (j=0;j<t;j++){
            long k;
            vector<int> x;
            if ((k=intersect(set1[j],c[i],dset1[j],d[i],x))<q) continue;

            if (k==d[i]) ok=0;
            if (k==dset1[j]){
                set2[j].push_back(i);
                dset2[j]++;
            }else{
//                assert(dd<max2);                
                long l;
                for (l=0;l<k;l++) set1[dd].push_back(x[l]);dset1[dd]=k;
                for (l=0;l<dset2[j];l++) set2[dd].push_back(set2[j][l]);
                
                dset2[dd]=dset2[j]+1;
                set2[dd].push_back(i);
                del[dd]=1;del[j]=0;
                dd++;                
            }
            ok2=0;
        }
        if (ok){
                assert(dd<max2);
            
                long l,k=d[i];
                for (l=0;l<k;l++) set1[dd].push_back(c[i][l]);dset1[dd]=k;
                dset2[dd]=1;
                set2[dd].push_back(i);
                if (ok2) del[dd]=1;else del[dd]=0;
                dd++;
        }
    }        
}

// is A a subset of B ?
long subset(vector<int> a,long m,vector<int> b,long n){
    int i,j=0;
    for (i=0;i<m;i++){
        while(j<n&&b[j]>a[i]) j++;
        if (j==n||b[j]<a[i]) return -1;
    }
    if (j==n-1&&m==n) return 0;
    return 1;
}

long check(long x){
    long i;
    for (i=0;i<dd;i++){
        long j=subset(set2[x],dset2[x],set2[i],dset2[i]);
        if (j==1) return 0;
        if (j==0&&i<x) return 0;
    }
    return 1;
}

char seq[max];
subpos c[max*2];
subpos d;

void print(){
    long i,j;
    for (i=0;i<dd;i++)
        if (del[i]&&dset2[i]>1&&check(i))
        {
            char s[max];
            int n=strlen(seq);
            for (j=0;j<n;j++) s[j]='N';
            
            for (j=0;j<dset2[i];j++)
                for (int l=0;l<moflen;l++) s[set2[i][j]+l]=seq[set2[i][j]+l];

            int dd=0;
            for (j=0;j<n;j++)
                if (s[j]!='N') dd++;                        
            s[n]=0;
            if (dd>=cov){         
                int l=0;
                int k=strlen(s)-1;
                while(l<k&&s[l]=='N') l++;
                while(l<k&&s[k]=='N') k--;
                
                for (int i1=l;i1<=k;i1++) cout<<s[i1];
                cout<<endl;
                cout<<dset1[i]+1<<endl;
                cout<<d.seq<<" "<<d.pos+l<<endl;
                for (j=0;j<dset1[i];j++) cout<<c[set1[i][j]].seq<<" "<<c[set1[i][j]].pos+l<<endl;
            }
        }
}

long a[max*2][max3];

void mining(){
    int q;
    cin>>q;
    long b[max*2];
    int n,j;
    
    long mx=0;
    int skip=0;
    dd=0;
    
    while(1){        
        cin>>seq;    
        if (strncmp(seq,"end",3)==0) break;
        
        cin>>moflen>>cov>>d.seq>>d.pos;
        n=0;
        while(1){
            assert(n<max);
            cin>>c[n].seq>>c[n].pos>>b[n];
            if (b[n]==-1) break;
            for (j=0;j<b[n];j++) cin>>a[n][j];            
            n++;            
        }
        
/*        if (n>=max){
            skip++;continue;
        }*/
        
        if (n>mx) mx=n;
        
        merge(a,b,n,q);
        print();
    }
    cout<<"end"<<endl;
    cout<<mx<<endl;
    cout<<"skip :"<<skip<<endl;
}

int main(){
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    time_t seconds,end;
    seconds = time (NULL);

    mining();

    end = time (NULL);
    cout<<"Total time: "<<end-seconds<<endl;  
}
