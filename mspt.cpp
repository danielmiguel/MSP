
#include <iostream>
#include<vector>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <string>
#include<time.h>



using namespace std;





#pragma region "Graph classes"

template<class T> struct Edge
{
    int v1;
    int v2;
    T cost;

    Edge<T>():v1(-1),v2(-1),cost(numeric_limits<T>::max()){}

    Edge<T>(int vert1, int vert2, T edgecost):v1(vert1),v2(vert2), cost(edgecost){}
 


	friend ostream & operator<<(ostream & str, const Edge<T> & val) { 
	  
		return str << "v1: "<< val.v1 << "; v2: "<< val.v2 << "; cost: "<< val.cost << "; ##";
	}

};

template<class T>
inline bool operator==(const Edge<T> &left, const Edge<T> &right)
{
	return (left.v1 == right.v1 && left.v2 == right.v2) || (left.v2 == right.v1 && left.v1 == right.v2);    
}

template<class T>
inline bool operator<(const Edge<T> &left, const Edge<T> &right)
{	
	return left.cost < right.cost;
}

template<class T>
inline bool operator>(const Edge<T> &left, const Edge<T> &right)
{	
	return left.cost > right.cost;
}

template<class T> struct Path
{
    vector<int> vertices;  
    T cost;
    bool exists;
    int vertex1,vertex2;

    Path(bool pathexists, vector<int> vertexset, T pathCost,int pathVertex1, int pathVertex2):exists(pathexists), vertices(vertexset), cost(pathCost), vertex1(pathVertex1), vertex2(pathVertex2){}

    Path():exists(false), cost(numeric_limits<T>::max()), vertex1(-1), vertex2(-1){ }
};


template<class T> class Graph
{

    public:
        Graph<T>(int verticesamount):verticesAmount(verticesamount)
        {          
            this->data = new vector<vector<Edge<T> > >(verticesamount);
        }
		Graph<T>(string datapath)
        {          
			cout<<"creating graph"<<endl;
			if(!this->LoadGraphFromFile(datapath))
				throw invalid_argument("The file does not exist or is not in correct format");
			cout<<"graph created"<<endl;
        }

        ~Graph<T>()
        {
            delete this->data;
        }
        inline int VerticesCount()
        {
            return this->verticesAmount;
        }

        inline int EdgesCount()
        {
            int sum = 0;
            for(int i = 0; i < data->size(); i++)
                sum += data->at(i).size();   

            return sum/2;
        }

        void AddEdge(int vertex1, int vertex2, T weight)
        {
            Edge<T> e = Edge<T>(vertex1, vertex2, weight);

            typename vector<Edge<T> >::iterator it;

            it = find( this->data->at(vertex1).begin(), this->data->at(vertex1).end(), e);

            if(it == this->data->at(vertex1).end())
                this->data->at(vertex1).push_back(e);

            it = find( this->data->at(vertex2).begin(), this->data->at(vertex2).end(), e);

            if(it == this->data->at(vertex2).end())
                this->data->at(vertex2).push_back(Edge<T>(vertex2, vertex1, weight));

        }

        void AddEdge(Edge<T> e)
        {
            this->AddEdge(e.v1, e.v2, e.cost);
        }

        void DeleteEdge(int vertex1, int vertex2)
        {
            Edge<T> e = Edge<T>(vertex1, vertex2, T());

            typename vector<Edge<T> >::iterator it;

            it = find( this->data->at(vertex1).begin(), this->data->at(vertex1).end(), e);

            if(it != this->data->at(vertex1).end()){
                this->data->at(vertex1).erase(it,it+1);           
                vector<Edge<T> >(this->data->at(vertex1)).swap(data->at(vertex1));
            }

            it = find( this->data->at(vertex2).begin(), this->data->at(vertex2).end(), e);

            if(it != this->data->at(vertex2).end()){
                this->data->at(vertex2).erase(it,it+1);
                vector<Edge<T> >(this->data->at(vertex2)).swap(data->at(vertex2));
            }
        }

        vector<int> Neighbors(int vertex)
        {
            vector<int> result;
            for(int i = 0; i < this->data->at(vertex).size(); i++)
                result.push_back(this->data->at(vertex)[i].v2);

            return result;
        }

        vector<Edge<T> > NeighborsEdges(int vertex)
        {
            if(vertex >= this->data->size())
                return vector<Edge<T> >(0);

            return this->data->at(vertex);
        }

        bool ContainsEdge(Edge<T> e)
        {
            typename vector<Edge<T> >::iterator it;
            it = find( this->data->at(e.v1).begin(), this->data->at(e.v1).end(), e);
            if(it != this->data->at(e.v1).end())
                return true;

            it = find( this->data->at(e.v2).begin(), this->data->at(e.v2).end(), e);
            if(it != this->data->at(e.v2).end())
                return true;
            	
			return false;
        }

		T EdgeCost(int vertex1, int vertex2)
		{
			typename vector<Edge<T> >::iterator it;
			it = find( this->data->at(vertex1).begin(), this->data->at(vertex1).end(), Edge<T>(vertex1, vertex2, T()));
			if(it != this->data->at(vertex1).end())
				return it->cost;                

            it = find( this->data->at(vertex2).begin(), this->data->at(vertex2).end(), Edge<T>(vertex1, vertex2, T()));
            if(it != this->data->at(vertex2).end())
               return it->cost; 

			return numeric_limits<T>::max();
		}

		bool ContainsEdge(int vertex1, int vertex2)
        {
			return ContainsEdge(Edge<T>(vertex1, vertex2, T()));
		}

		vector<Edge<T> > GetEdges()
		{
			vector<Edge<T> > result;

			Edge<T> current;

			for(int i = 0; i < this->data->size(); i++)
			{
				for(int j = 0; j <  this->data->at(i).size(); j++)
				{
					current = this->data->at(i).at(j);

					if(current.v2 > i)
						result.push_back(current);
				}
			}

			return result;
		}

		/*
		*Gets the minimum spanning tree using Kruskal algorithm
		*sum:Pointer to the sum of the edges of the tree(Minimum Spanning Tree)
		*/
		vector<Edge<T> > MSPT_KRUSKAL(T *sum)
		{
			*sum = 0;
			vector<Edge<T> > result;


			//vector to represent the forest.
			//At the beginning, each index represents a tree of only one node
			vector<int> forestIndex(this->verticesAmount, -1);

			
			cout<<"starts getting edges"<<endl;

			//edges: vector with the edges of the tree
			vector<Edge<T> > edges = this->GetEdges();	

			cout<<"starts sorting the edges list"<<endl;			

			//sorting the edges using the cost(operator < and > overloaded)			
			sort(edges.begin(), edges.end());

			cout<<"edges sorted"<<endl;


			int edgecounts = 0;

			Edge<T> current; 
      
			
			int v1 = 0;//vertex 1
			int v2 = 0;//vertex 2
			int root1 = -1;// root of the tree that vertex1 belongs to
			int root2 = -1;//root of the tree that vertex2 belongs to

			int edgessize = edges.size();
			int forestsize = forestIndex.size();
			
			int currentIndex = 0;
			int dp1 = 0;
			int dp2 = 0;

			//kruskal algorithm

			//stop condition: we iterated over all the edges or we have verticesAmount - 1 edges
			
			for(int i =0; i < edgessize && edgecounts < this->verticesAmount - 1; i++)
			{
				current = edges[i];
				v1 = current.v1;//vertex1
				v2 = current.v2;//vertex2

				root1 = -1;
				root1 = -1;
				 dp1 = 0;
				 dp2 = 0;

				 //calculating the root of the vertex1 by getting the parent until we reach the node that does not have parent
				currentIndex = v1;
				while(forestIndex[currentIndex] != -1)
				{
					currentIndex = forestIndex[currentIndex];
					dp1++;
				}
				root1 = currentIndex;


				//second loop to try to make all the parents of v1 as roo1 in order to make the tree more wide and less deep
				currentIndex = v1;
				while(forestIndex[currentIndex] != -1)
				{
					currentIndex = forestIndex[currentIndex];					
					if(forestIndex[currentIndex] != -1)
						forestIndex[currentIndex] = root1;
				}
				
				
				
				//calculating the root of the vertex2 by getting the parent until we reach the node that does not have parent
				currentIndex = v2;
				while(forestIndex[currentIndex] != -1)
				{
					currentIndex = forestIndex[currentIndex];
					dp2++;
				}
				root2 = currentIndex;	


				//second loop to try to make all the parents of v2 as roo2 in order to make the tree more wide and less deep
				currentIndex = v2;
				while(forestIndex[currentIndex] != -1)
				{
					currentIndex = forestIndex[currentIndex];
					if(forestIndex[currentIndex] != -1)
						forestIndex[currentIndex] = root2;
				}

				//if both vertices have the same root, we don't need to do nothing because they are already in the same forest
				
				if(root1 != root2)//if the vertices are in different forests we merge this forest by setting the tree with the smaller depth as child of the tree
				//with the greatest depth. By doing this we make sure that we keep the trees with the shortest depth as possible. If we don't do this, calculating the roots of
				//vertex1 and vertex2 would be very time-expensive.
				{		
					if(dp1 < dp2)
					{
						forestIndex[root1] = root2;
					}
					else
					{
						forestIndex[root2] = root1;
					}					
					
					result.push_back(current);
					*sum += current.cost;
					edgecounts++;
				}
		
			}
				
			

			return result;
		}

    private:           
        int verticesAmount;

        vector<vector<Edge<T> > > *data;
		/*
		*Create graph edges from a file.
		*Format:
		*first line: vertices amount
		*rest of the lines: vertex1 vertex2 cost
		*/
	
		
		bool LoadGraphFromFile(string path)
		{
			
			ifstream infile;
			infile.open(path.c_str());
			std::string line;
			
			//getting the vertices amount. First line
			getline(infile, line);

			istringstream iss(line);
			
			int verticesAmount = 0;
			if (!(iss >> verticesAmount )) { infile.close();return false; } // error
			
			//now in verticesAmount we have the vertices amount

			//creating the edges vector
			this->data = new vector<vector<Edge<T> > >(verticesAmount);
			this->verticesAmount = verticesAmount;

			int vertex1, vertex2;
			T cost;
						
			string::size_type sz;
			string::size_type sz2;

			//for each line we create the edges
			while (getline(infile, line))
			{

				//obtaining the values of vertex1, vertex2 and cost from the line
				
				sscanf(line.c_str(), "%d %d %lf", &vertex1, &vertex2, &cost);
				

				//Adding the edge. For speed purpose we add the edge directly in the vector. 
				//If we want more checking, we can use this->AddEgde
				//this->AddEdge(vertex1, vertex2, cost);
				this->data->at(vertex1).push_back(Edge<T>(vertex1, vertex2, cost));				
			}

			infile.close();
			return true;
		}

		

};

#pragma endregion


int main(int argc, char *argv[])
{
	try{
		
		 double seconds;
		 double cost = 0;
		
		Graph<double> *graph = new Graph<double>("F:\\graphsample6.txt");				
	
		
		clock_t cl = clock();

		vector<Edge<double> > mspt = graph->MSPT_KRUSKAL(&cost);


		/*for each (Edge<double>  var in mspt)
		{
			cout<<var<<endl;
		}*/
		

		cout<<endl;
		cout<<"cost is: "<<cost<<endl;
		cout<<"It took "<<(static_cast<double>(clock() - cl)/1000)<<" seconds."<<endl;

	}
    catch(exception& e)
    {
        cout << e.what() << endl;
    }
	cin.get();
	return 0;
}
