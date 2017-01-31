#include <cmath>
#include <algorithm>
#include <ctime>
#include <iostream>
#include <string.h>
#include <getopt.h>
#include <sstream>
#include <vector>
#include <queue>
#include <cstdio>
#include <cstdlib>
#include <iomanip>

using namespace std;

struct node{
  int id;
  node *p;
  bool k;
  int x;
  int y;
};

struct edge{
  node *head;
  node *tail;
  double d;
};

struct TSPcomp{
  bool operator()(const edge* x, 
			  const edge* y) const{
	if(y->d == -1) return 1;
	else if(x->d == -1) return 0;
	else return (x->d < y->d);
  }
};

struct TSP{
	node *N;
	TSP* p;
	bool k;
	int level;
	double d;
};

bool isout(const node* n, 
			  const int &lx, const int &ly, 
			  const int &rx ,const int &ry) {
	return (n->x <= lx) || (n->x >= rx) 
			|| (n->y <= ly) || (n->y >= ry);
}

bool isin(const node* n, 
			  const int &lx, const int &ly, 
			  const int &rx ,const int &ry) {
	return (n->x >= lx) && (n->x <= rx) 
			&& (n->y >= ly) && (n->y <= ry);
}

void ope(TSP* Te, vector<TSP*>& result, vector<TSP*>& T, vector<vector<edge*>>& L, double& Tweight, double d);

int main(int argc, char* argv[])
{
  cout << std::setprecision(2);
  cout << std::fixed;
  string line;
  string ex_type;
  extern char *optarg;
  int arg;
  const char* const short_options = "m:h";
  const struct option long_options[] = {
	{"mode",	1,	0,	'm'},
	{"help",	0,	0,	'h'},
	{0,			0,	0,	 0 }
  };

  do{
	arg = getopt_long(argc, argv, short_options, long_options, NULL);
	switch (arg){
	case 'm':{
	  ex_type = optarg;
	  break;
	}
	case 'h':
	  exit(0);
	  break;
	case -1:
	  break;
	default:
	  abort();
	}
  }while(arg != -1);

////////////////////////////////////////////////////////////////////////
  int lx = 0, ly = 0, rx = 0, ry = 0;
  int PNum;
  istringstream iss;
  if(ex_type == "MST"){
    getline(cin, line);
    iss.str(line);
	iss >> lx >> ly >> rx >> ry;
	iss.clear();
    getline(cin, line);
    iss.str(line);
    iss >> PNum;
	iss.clear();
    vector<node*> vMST;
	vector<double> key;
	vector<node*> Parent;
    int count = 0;
	while(getline(cin, line) && count < PNum){
	  iss.str(line);
	  node* n = new node;
      n->id = count;
      n->k = 0;
	  iss >> n->x >> n->y;
	  vMST.push_back(n);
	  key.push_back(-1);
	  Parent.push_back(nullptr);
	  iss.clear();
	  count++;
    }
    int index = 0;
	int EN = 0;
    double Tweight = 0;
	int min = -1;
    node* minP;
	if(PNum != 0) vMST[index]->k = 1;
    while(EN != PNum && PNum != 0){
	  if(min == -1 && EN != 0) break;
	  if(min != -1){
		vMST[min]->k = 1;
		vMST[min]->p = minP;
		Tweight += key[min];
		index = min;
	  }
	  min = -1;
	  for(int j = 0; j < PNum; j++){
		if(j == index) continue;
		if(vMST[j]->k) continue;
		if((min == -1 && key[j] != -1) || 
		   (min != -1 && key[j] < key[min] && key[j] != -1)){
		  min = j;
		  minP = Parent[j];
		}
		if((isout(vMST[index], lx, ly, rx, ry) && 
		   isout(vMST[j], lx, ly, rx, ry)) ||
		   (isin(vMST[index], lx, ly, rx, ry) && 
		   isin(vMST[j], lx, ly, rx, ry))){
		  double newd = 0;
		  newd = sqrt(pow((vMST[index]->x - vMST[j]->x), 2.0) + 
					 pow((vMST[index]->y - vMST[j]->y), 2.0));
		  if((min != -1 && key[min] > newd) ||
			 (min == -1)){
			min = j;
			minP = vMST[index];
		  }
		  if(newd < key[j] || key[j] == -1){
		    key[j] = newd;
			Parent[j] = vMST[index];
		  }
		}
      }
      EN++;
    }
	if(EN < PNum || PNum == 0) {
	  fprintf(stderr, "Cannot construct MST!\n");
	  exit(1);
	}
	cout << Tweight << "\n";
	for(unsigned int i = 1; i < vMST.size(); i++){
	  int a = vMST[i]->id < vMST[i]->p->id ? 
			  vMST[i]->id : vMST[i]->p->id;
	  int b = vMST[i]->id < vMST[i]->p->id ? 
			  vMST[i]->p->id : vMST[i]->id;
	  cout << a << " " << b << "\n";
	}
  }
  if(ex_type == "OPTTSP"){
    getline(cin, line);
    getline(cin, line);
    iss.str(line);
    iss >> PNum;
	iss.clear();
    vector<node*> Ver;
    vector<TSP*> T;
	vector<TSP*> result;
	vector<vector<edge*>> L;
	double Tweight = -1;
    int count = 0;
	while(getline(cin, line) && count < PNum){
	  iss.str(line);
	  node* n = new node;
      n->id = count;
      n->k = 0;
	  iss >> n->x >> n->y;
	  Ver.push_back(n);
	  vector<edge*> newTv;
	  L.push_back(newTv);
	  iss.clear();
	  count++;
    }
	for(int i = 0; i < PNum; i++) {
	  for(int j = i + 1; j < PNum; j++) {
		double newd = 0;
		newd = sqrt(pow((Ver[i]->x - Ver[j]->x), 2.0) + 
			   pow((Ver[i]->y - Ver[j]->y), 2.0));
		edge* newedge1 = new edge;
		edge* newedge2 = new edge;
		newedge1->head = Ver[i];
		newedge1->tail = Ver[j];
		newedge1->d = newd;
		newedge2->head = Ver[j];
		newedge2->tail = Ver[i];
		newedge2->d = newd;
		L[i].push_back(newedge1);
		L[j].push_back(newedge2);
	  }
	  std::sort(L[i].begin(), L[i].end(), TSPcomp());
	  TSP* newTSP = new TSP;
	  newTSP->N = Ver[i];
	  newTSP->p = nullptr;
	  newTSP->k = 0;
	  newTSP->level = -1;
	  newTSP->d = 0;
	  T.push_back(newTSP);
	  result.push_back(newTSP);
	}
	T[0]->level = -1;
	T[0]->p = T[0];
	if(PNum > 1) {
	  ope(T[0], result, T, L, Tweight, 0);
	  cout << Tweight <<"\n";
	  cout << result[0]->N->id;
	  for(unsigned int i = 1; i < result.size(); i++)
		cout << " " << result[i]->N->id;
	  cout << "\n";
	}
  }
  if(ex_type == "FASTTSP"){
	clock_t curtime = clock();
    getline(cin, line);
    getline(cin, line);
    iss.str(line);
    iss >> PNum;
	iss.clear();
    vector<node*> Ver;
    vector<TSP*> T;
	vector<node*> result;
	double Tweight = 0;
    int count = 0;
	while(getline(cin, line) && count < PNum){
	  iss.str(line);
	  node* n = new node;
      n->id = count;
      n->k = 0;
	  iss >> n->x >> n->y;
	  Ver.push_back(n);
	  iss.clear();
	  count++;
    }
	double min = -1;
	int curnode = 0;
	int num = 1;
	int next = -1;
	result.push_back(Ver[0]);
	Ver[0]->k = 1;
	while(num < PNum){
	  min = -1;
	  next = -1;
	  for(int i = 0; i < PNum; i++){
		if(Ver[i]->k || i == curnode) continue;
		if(min == -1){
		  min = sqrt(pow((Ver[i]->x - Ver[curnode]->x), 2.0) + 
		   		pow((Ver[i]->y - Ver[curnode]->y), 2.0));
		  next = i;
		}
		else{
		  double newd = 0;
		  newd = sqrt(pow((Ver[i]->x - Ver[curnode]->x), 2.0) + 
		   		 pow((Ver[i]->y - Ver[curnode]->y), 2.0));
		  if(newd < min){
			min = newd;
			next = i;
		  }
		}
	  }
	  result.push_back(Ver[next]);
	  Ver[next]->k = 1;
	  curnode = next;
	  Tweight += min;
	  num++;
	}
	Tweight += sqrt(pow((Ver[0]->x - Ver[curnode]->x), 2.0) + 
		   	   pow((Ver[0]->y - Ver[curnode]->y), 2.0));
	num = 0;
	while(num < PNum - 1){
	  if((clock() - curtime) / CLOCKS_PER_SEC > 13) break;
	  int change = 0;
	  double diff = -1;
	  double old_d = 0, new_d = 0;
	  for(int i = num + 1; i < PNum; i++){
		int Aright = (i + 1) % PNum;
		int Bleft = (num == 0) ? PNum - 1 : num - 1;
		int Bright = num + 1;
		int Aleft = i - 1;
		if(i == num + 1 || (num == 0 && i == PNum - 1)){
		  old_d = sqrt(pow((result[i]->x - result[Aright]->x), 2.0) + 
		   		  pow((result[i]->y - result[Aright]->y), 2.0)) +
				  sqrt(pow((result[i]->x - result[Aleft]->x), 2.0) + 
		   		  pow((result[i]->y - result[Aleft]->y), 2.0)) +
				  sqrt(pow((result[num]->x - result[Bright]->x), 2.0) + 
		   		  pow((result[num]->y - result[Bright]->y), 2.0)) +
				  sqrt(pow((result[num]->x - result[Bleft]->x), 2.0) + 
		   		  pow((result[num]->y - result[Bleft]->y), 2.0)) -
				  sqrt(pow((result[num]->x - result[i]->x), 2.0) + 
		   		  pow((result[num]->y - result[i]->y), 2.0)) -
				  sqrt(pow((result[num]->x - result[i]->x), 2.0) + 
		   		  pow((result[num]->y - result[i]->y), 2.0));
		  new_d = sqrt(pow((result[i]->x - result[Bright]->x), 2.0) + 
		   		  pow((result[i]->y - result[Bright]->y), 2.0)) +
				  sqrt(pow((result[i]->x - result[Bleft]->x), 2.0) + 
		   		  pow((result[i]->y - result[Bleft]->y), 2.0)) +
				  sqrt(pow((result[num]->x - result[Aright]->x), 2.0) + 
		   		  pow((result[num]->y - result[Aright]->y), 2.0)) +
				  sqrt(pow((result[num]->x - result[Aleft]->x), 2.0) + 
		   		  pow((result[num]->y - result[Aleft]->y), 2.0));
		}
		else{
		  old_d = sqrt(pow((result[i]->x - result[Aright]->x), 2.0) + 
		   		  pow((result[i]->y - result[Aright]->y), 2.0)) +
				  sqrt(pow((result[i]->x - result[Aleft]->x), 2.0) + 
		   		  pow((result[i]->y - result[Aleft]->y), 2.0)) +
				  sqrt(pow((result[num]->x - result[Bright]->x), 2.0) + 
		   		  pow((result[num]->y - result[Bright]->y), 2.0)) +
				  sqrt(pow((result[num]->x - result[Bleft]->x), 2.0) + 
		   		  pow((result[num]->y - result[Bleft]->y), 2.0));
		  new_d = sqrt(pow((result[i]->x - result[Bright]->x), 2.0) + 
		   		  pow((result[i]->y - result[Bright]->y), 2.0)) +
				  sqrt(pow((result[i]->x - result[Bleft]->x), 2.0) + 
		   		  pow((result[i]->y - result[Bleft]->y), 2.0)) +
				  sqrt(pow((result[num]->x - result[Aright]->x), 2.0) + 
		   		  pow((result[num]->y - result[Aright]->y), 2.0)) +
				  sqrt(pow((result[num]->x - result[Aleft]->x), 2.0) + 
		   		  pow((result[num]->y - result[Aleft]->y), 2.0));
	    }
		if((diff == -1 && old_d > new_d) || (diff != -1 && diff < old_d - new_d)){
		  diff = old_d - new_d;
		  change = i;
		}
		//cout<<new_d<<" "<<old_d<<" A "<<result[i]->id<<" Aleft "<<result[Aleft]->id<<" Aright "<<result[Aright]->id<<endl;
		//cout<<diff<<" "<<" B "<<result[num]->id<<" Bleft "<<result[Bleft]->id<<" Bright "<<result[Bright]->id<<endl;
	  }
	  if(diff != -1){
		node* temp = result[change];
		result[change] = result[num];
		result[num] = temp;
		//cout << diff << endl;
		Tweight = Tweight - diff;
	  }
	  num++;
	}
	if(PNum > 1) {
	  cout << Tweight <<"\n";
	  int h = 0;
	  for(h = 0; h < (int) result.size(); h++){
		if(result[h]->id == 0) break;
	  }
	  cout << result[h]->id;
	  for(unsigned int i = h + 1; i < result.size() + h; i++)
		cout << " " << result[i % result.size()]->id;
	  cout << "\n";
	}
  }
  return 0;
}

void ope(TSP* Te, vector<TSP*>& result, vector<TSP*>& T, vector<vector<edge*>>& L, double& Tweight, double d){
  Te->level = Te->p->level + 1;
  Te->d = Te->p->d + d;
  Te->k = 1;
  //cout << "Te: " << "id: " <<Te->N->id << "level: "<<Te->level << "d: " <<Te->d <<"d: "<<d<<"k: "<<Te->k << endl;
  if(Te->level == (int) T.size()-1) {
	//update result
	result[0] = T[0];
	TSP* track = Te;
	for(int i = 1; i < (int) T.size(); i++){
	  //cout<<"track: "<<track->N->id<<endl;
	  result[i] = track;
	  track = track->p;
	}
	double newd = 0;
	newd = sqrt(pow((Te->N->x - T[0]->N->x), 2.0) + 
		   pow((Te->N->y - T[0]->N->y), 2.0));
	Tweight = newd + Te->d;
	//cout<<"first if: "<<Tweight<<" "<<newd<<"Te-d"<<Te->d<<endl;
	Te->k = 0;
	return;
  }
  for(int i = 0; i < (int) T.size() - 1; i++) {
	//cout<<"i="<<i<<endl;
	if(!T[L[Te->N->id][i]->tail->id]->k){
	  double lowerbound = 0;
	  for(int j = 0; j < (int) T.size(); j++) {
	    //cout<<"\t"<<"j="<<j<<endl;
	    if(!T[j]->k || j == 0) {
		  int l = 0;
		  while(T[L[j][l]->tail->id]->k && 
			    L[j][l]->tail->id != 0) l++;
		  //cout<<"\tl after ++1="<<l<<"head="<<L[j][l]->head->id<<"tail="<<L[j][l]->tail->id<<endl;
		  lowerbound += L[j][l]->d / 2;
		  l++;
		  if((j != L[Te->N->id][i]->tail->id) && j != 0) {
		    while(T[L[j][l]->tail->id]->k && 
				  L[j][l]->tail->id != 0) l++;
		    //cout<<"\tl after ++2="<<l<<"head="<<L[j][l]->head->id<<"tail="<<L[j][l]->tail->id<<endl;
		    lowerbound += L[j][l]->d / 2;
		  }
	    }
	  }
	  //cout<<"L[Te->N->id][i]="<<L[Te->N->id][i]->head->id<<"tail="<<L[Te->N->id][i]->tail->id<<endl;
	  //cost = (cost + L[Te->N->id][i]->d) / 2 + Te->d;
	  double cost = lowerbound + L[Te->N->id][i]->d + Te->d;
	  //cout<<"cost="<<cost<<"Tweight="<<Tweight<<endl;
	  //if((cost >= Tweight) && (Tweight != -1)) {cout<<"break"<<endl;break;}
	  if((cost >= Tweight) && (Tweight != -1)) break;
	  T[L[Te->N->id][i]->tail->id]->p = Te;
	  ope(T[L[Te->N->id][i]->tail->id], result, T, L, Tweight, L[Te->N->id][i]->d);
	}
  }
  Te->k = 0;
  return;
}
