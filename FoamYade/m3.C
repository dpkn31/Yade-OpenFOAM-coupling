
/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    meshAndField

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "meshSearch.H"
#include <vector>
#include <algorithm> 

class meshpt {
  
  public:
    meshpt(const double* x, const double* y, const double* z, int cid){
    pt.push_back(x); pt.push_back(y);  pt.push_back(z); id = cid;  } 
    std::vector<const double*> pt; 
    int id; 
   virtual ~meshpt(){}; 
}; 




class kdNode {
  
  public:  
    kdNode(meshpt point) :p(point) {left=NULL; right=NULL; }
    meshpt p; 
    kdNode* left; 
    kdNode* right;

    virtual ~kdNode(){};
}; 


class querypt {


  public: 

    querypt(meshpt p1) : p(p1) {}
    meshpt p;  
    double best_dist=1e9; 
    std::vector<int> listids; 
    virtual ~querypt(){} 
    kdNode* bestnode; 
}; 




class kdTree { 
 public:
  kdTree(){}; 
  kdNode* root; 
  void build_tree(std::vector<meshpt>&); 
  kdNode* recursive_build_tree(std::vector<meshpt>&, int); 
  int nearestCell(const vector& );  
  const int ndim=3; 
  typedef std::vector<meshpt>::size_type vec_sz;
  void get_median(std::vector<meshpt>& ,const int&);
  kdNode* recursive_nearest_cell(kdNode* , const meshpt& , kdNode* , double& , int);
  double distance(const meshpt& , const meshpt& );  
  int numlevels; 
  virtual ~kdTree() {};
  void traversTree(); 
  void _traversTree(kdNode*); 
   
}; 



class cmpvec{

  public: 
    cmpvec(int axis) : a(axis) {}; 

    bool operator () (const meshpt& p1, const meshpt& p2)
    {
      return *(p1.pt[a]) < *(p2.pt[a]); 

    }

    int a; 
}; 




void kdTree::build_tree(std::vector<meshpt>& ptlist) 
{
  numlevels =0; 
  root = recursive_build_tree(ptlist, 0);
   

}


kdNode* kdTree::recursive_build_tree(std::vector<meshpt>& ptlist, int depth) 
{

  if (! ptlist.size()){return NULL; } 
  numlevels = numlevels+1; 
  int axis = depth%ndim; 

 
  get_median(ptlist, axis);
  vec_sz md = ptlist.size()/2;  
  kdNode* node = new kdNode(ptlist[md]);

  std::vector<meshpt> pv1 = std::vector<meshpt>(ptlist.begin(), ptlist.begin()+md); 
  std::vector<meshpt> pv2 = std::vector<meshpt>(ptlist.begin()+md+1, ptlist.end()); 


  node->left = recursive_build_tree(pv1, depth+1);
  node->right = recursive_build_tree(pv2, depth+1);  

  return node;  


}

void kdTree::traversTree()
{
  _traversTree(root); 

}


void kdTree::get_median(std::vector<meshpt>& ptlist, const int& axis) 
{
//  std::sort(ptlist.begin(), ptlist.end(), cmpvec(axis)); 
   vec_sz md = ptlist.size()/2;   
  std::nth_element(ptlist.begin(), ptlist.begin()+md, ptlist.end(), cmpvec(axis)); 

}


double kdTree::distance(const meshpt& p2, const meshpt& p1)
{
  double dist = 0.0; 
  for (unsigned int i=0; i!= p1.pt.size();++i)
  {
    double ds = *(p1.pt[i])-*(p2.pt[i]); 
    dist += ds*ds; 
  }
  return dist; 

}

int kdTree::nearestCell(const vector& p)
{
  
//  kdNode* bnode;  
  meshpt v(& p.x(), &p.y(), &p.z(), -1);
     
  double dist = distance(root->p,  v); 
  kdNode* bnode =recursive_nearest_cell(root, v, root,dist,0); 
  dist =0.0; 
  dist = distance(bnode->p, v); 
  dist = Foam::sqrt(dist); 
  std::cout << "best distance = " << dist << std::endl; 
  return bnode->p.id; 
}

//
kdNode* kdTree::recursive_nearest_cell(kdNode* node, const meshpt& v, kdNode* best, double& best_dist,int depth)
{


 
  if (node ==NULL) return NULL;  
//  std::cout << "nearest pt = " << *(node ->p.pt[0]) <<  "   " << *(node ->p.pt[1]) << "  " << *(node ->p.pt[2]) << std::endl;  
//  std::cout << "best dist = " << qpoint.best_dist << std::endl;
//

  kdNode* best1=best;  
  double dist_l = best_dist; 
  double distsq = distance(node->p, v); 

   if (distsq < best_dist){ 
      dist_l = distsq;
      best1 = node; 

   }

   int axis = depth%ndim; 
   double df = *(node ->p.pt[axis])-*(v.pt[axis]);
   double df2 = df*df; 


   kdNode* next; 
   kdNode* other; 

//   bool 

   if (df > 0.0){
     next = node->left; 
     other =node->right;

      
   } else { 
     next = node-> right; 
     other = node-> left;
   }

  kdNode* nextN = recursive_nearest_cell(next, v,best1, dist_l, depth+1); 
  
  if (nextN != NULL){
    distsq =distance(nextN->p, v); 
    if (distsq < dist_l){ 
      dist_l = distsq;
      best1 = nextN;  
  }
}

  if(df2 < dist_l ){  

  kdNode* nextN = recursive_nearest_cell(other, v,best1, dist_l, depth+1); 
  if (nextN != NULL){

    distsq =distance(nextN->p, v); 
    if (distsq < dist_l){ 
      dist_l = distsq;
      best1 = nextN;  
    }
  }

}
return best1;








}


void kdTree::_traversTree(kdNode* node) 
{
 
  if(node==NULL){return; }
  _traversTree(node->left); 

  _traversTree(node ->right);   


}


//lass kdNode { 
//
//  const std::vector<double*> 
//
//}; 
//



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"


  std::vector<meshpt>  cellC;  
  forAll(mesh.C(), cellI)
//meshpt c_p(&mesh.C()[id].x(), &mesh.C()[id].y(), &mesh.C()[id].z(), id); 
    cellC.push_back(meshpt(&mesh.C()[cellI].x(), &mesh.C()[cellI].y(), &mesh.C()[cellI].z(), cellI));

 kdTree aTree; 
 aTree.build_tree(cellC);

 vector v(0.5, 0.5, 0.5); 
 int id = aTree.nearestCell(v); 
 label pid = mesh.findCell(v); 

 Info << "point from my shit = " << mesh.C()[id].x() << "  " << mesh.C()[id].y() << "  " << mesh.C()[id].z() << "  id = "<< id  << endl; 
 Info << "point from foam shit = " << mesh.C()[pid].x() << "  " << mesh.C()[pid].y() << "  " << mesh.C()[pid].z() << "  id = "<< pid << endl; 


return 0 ; 

  

  
}


// ************************************************************************* //
