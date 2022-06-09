#include<iostream> 
#include <list> 
#include <stdio.h> 
#include <stdlib.h> 
#include <time.h> 
#include <algorithm>
#include <random>
#include <cmath>
#include <ctime>
#include <chrono>
#include <fstream>
#include <string>
#include <algorithm>
using namespace std;

default_random_engine generator;

class Graph
{
	int V;
	list<int> *adj;
public:
	Graph(int V);
	void addEdge(int v, int w);
	int printVertexCover();
	void populate(int E);
	int exactVertexCover();
	void returnVertexCover(vector<int> &);
	void printGraph();
	bool isVertexCover(vector<int> &);
};

void Graph::printGraph() {
	for (int i = 0; i < V; i++) {
		list<int>::iterator k;
		for (k = adj[i].begin(); k != adj[i].end(); ++k)
		{
			int v = *k;
			cout << i << "---" << v << endl;
		}
	}
}

Graph::Graph(int V)
{
	this->V = V;
	adj = new list<int>[V];
}

void Graph::addEdge(int v, int w)
{
	adj[v].push_back(w);
	adj[w].push_back(v);
}

void makeCombiUtil(vector<vector<int> >& ans,
	vector<int>& tmp, int n, int left, int k)
{

	if (k == 0) {
		ans.push_back(tmp);
		return;
	}


	for (int i = left; i <= n; ++i)
	{
		tmp.push_back(i);
		makeCombiUtil(ans, tmp, n, i + 1, k - 1);


		tmp.pop_back();
	}
}


vector<vector<int> > makeCombi(int n, int k)
{
	vector<vector<int> > ans;
	vector<int> tmp;
	makeCombiUtil(ans, tmp, n, 1, k);
	return ans;
}

int Graph::exactVertexCover() {
	vector<int> vertexCover;
	for (int i = 1; i <= V; i++) {
		vector<vector<int>> combinations = makeCombi(V, i);
		for (int k = 0; k < combinations.size(); k++) {
			for (int j = 0; j < combinations[k].size(); j++) {
				list<int>::iterator i;
				vertexCover.push_back(combinations[k][j]);
				for (i = adj[combinations[k][j] - 1].begin(); i != adj[combinations[k][j] - 1].end(); ++i) {
					int v = *i;
					if (find(vertexCover.begin(), vertexCover.end(), v) == vertexCover.end()) {
						vertexCover.push_back(v);
					}
				}
			}
			if (vertexCover.size() == V) {
				/*cout << "Vertex Cover: " << endl;
				for (int l = 0; l < combinations[k].size(); l++) {
					cout << combinations[k][l] - 1 << endl;
				}*/
				return combinations[k].size();
			}
			vertexCover.clear();
		}

	}
	return -1;
}

int Graph::printVertexCover()
{
	//int vertexCoverSize = 0;
	bool * visited = new bool[V];
	for (int i = 0; i < V; i++) {
		visited[i] = false;
	}

	vector<int> vertexCover;

	list<int>::iterator i;


	for (int u = 0; u < V; u++)
	{

		if (visited[u] == false)
		{

			for (i = adj[u].begin(); i != adj[u].end(); ++i)
			{
				int v = *i;
				if (visited[v] == false)
				{

					visited[v] = true;
					visited[u] = true;
					vertexCover.push_back(v);
					vertexCover.push_back(u);
					for (i = adj[u].begin(); i != adj[u].end(); ++i) {
						int l = *i;
						visited[l] = true;
					}

					for (i = adj[v].begin(); i != adj[v].end(); ++i) {
						int l = *i;
						visited[l] = true;
					}

					break;
				}
			}
		}
	}
	/*	for (int i = 0; i < V; i++)
		if (visited[i]) {
			cout << i << " ";
			vertexCoverSize++;
		}
	cout << "Vertex Cover Size: " << vertexCoverSize << endl;
	return vertexCoverSize;
	*/
	return vertexCover.size();
}

bool Graph::isVertexCover(vector<int> & vertexCover) {

	bool * visited = new bool[V];
	for (int i = 0; i < V; i++) {
		visited[i] = false;
	}

	for (int k = 0; k < vertexCover.size(); k++) {
		list<int>::iterator i;
		visited[vertexCover[k]] = true;
		for (i = adj[vertexCover[k]].begin(); i != adj[vertexCover[k]].end(); ++i) {
			int v = *i;
			visited[v] = true;
		}
	}

	for (int i = 0; i < V; i++) {
		if (visited[i] == false) {
			return false;
		}
	}
	return true;
}

void Graph::returnVertexCover(vector<int> & vertexCover) {
	bool * visited = new bool[V];
	for (int i = 0; i < V; i++) {
		visited[i] = false;
	}

	list<int>::iterator i;


	for (int u = 0; u < V; u++)
	{

		if (visited[u] == false)
		{

			for (i = adj[u].begin(); i != adj[u].end(); ++i)
			{
				int v = *i;
				if (visited[v] == false)
				{

					visited[v] = true;
					visited[u] = true;
					break;
				}
			}
		}
	}
	for (int i = 0; i < V; i++)
		if (visited[i]) {
			vertexCover.push_back(i);
		}
}

void Graph::populate(int E) {

	int fromV = 0;
	int toV = 0;
	uniform_int_distribution<int> distribution(0, V - 1);
	for (int i = 0; i < E; i++) {
		fromV = distribution(generator);
		toV = distribution(generator);

		while (fromV == toV || find(adj[toV].begin(), adj[toV].end(), fromV) != adj[toV].end()) {
			fromV = distribution(generator);
			toV = distribution(generator);
		}

		addEdge(fromV, toV);
	}

	list<int>::iterator i;
	for (int u = 0; u < V; u++)
	{
		if (adj[u].empty()) {
			toV = distribution(generator);
			addEdge(u, toV);
		}
	}
}

float calculateSD(vector<float> & data) {
	float sum = 0.0, mean, stadardDeviation = 0.0;
	for (int i = 0; i < data.size(); i++) {
		sum += data[i];
	}
	mean = sum / data.size();
	for (int j = 0; j < data.size(); j++) {
		stadardDeviation += pow(data[j] - mean, 2);
	}
	stadardDeviation = sqrt(stadardDeviation / data.size());
	return stadardDeviation;
}

float calculateStandardError(float standardDeviation, int N) {
	return standardDeviation / sqrt(N);
}

void getRunningTime(vector<float> & runningTimes, ofstream & myFile) {
	float totalTime = 0.0;
	int N = runningTimes.size();
	for (int i = 0; i < N; i++) {
		totalTime += runningTimes[i];
	}

	float standardDeviation = calculateSD(runningTimes);


	float m = totalTime / N;

	const float tval90 = 1.645;

	const float tval95 = 1.96;

	float sm = calculateStandardError(standardDeviation, N);

	float upperMean90 = m + tval90 * sm;
	float lowerMean90 = m - tval90 * sm;

	float upperMean95 = m + tval95 * sm;
	float lowerMean95 = m - tval95 * sm;

	/*cout << "Mean Time: " << m << " ms" << endl;
	cout << "SD: " << standardDeviation << endl;
	cout << "Standard Error: " << sm << endl;
	cout << "%90-CL: " << upperMean90 << " - " << lowerMean90 << endl;
	cout << "%95-CL: " << upperMean95 << " - " << lowerMean95 << endl;
	*/

	myFile << to_string(m) + ",";
	myFile << to_string(standardDeviation) + ",";
	myFile << to_string(sm) + ",";
	myFile << to_string(lowerMean90) + " - " + to_string(upperMean90) + ",";
	myFile << to_string(lowerMean95) + " - " + to_string(upperMean95) + ",\n";
	runningTimes.clear();
}

// Driver program to test methods of graph class 
int main()
{
	int option;

	cout << "Enter 1 for Running Time Analysis: " << endl;
	cout << "Enter 2 for Quality Analysis: " << endl;
	cout << "Enter 3 for Correctness Analysis:" << endl;
	cin >> option;
	if (option == 1) {
		ofstream myfile;
		int E, V, n;
		int N[] = { 100, 1000, 5000, 10000 };
		for (int i = 0; i < 4; i++) {
			n = N[i];
			myfile.open("EdgeFixed-" + to_string(n) + "-Iterations.csv");
			myfile << "Size(V),";
			myfile << "Mean Time(ms),";
			myfile << "Standard Deviation,";
			myfile << "Standard Error,";
			myfile << "% 90 - CL,";
			myfile << "% 95 - CL\n";

			V = 100;
			E = 100;

			for (int k = 0; k < 10; k++) {
				int new_vertex = V + (k * 100);
				myfile << to_string(new_vertex) + ",";
				vector<float> runningTimes;

				for (int i = 0; i < n; i++) {
					Graph G(new_vertex);
					G.populate(E);
					auto t_start = std::chrono::high_resolution_clock::now();
					G.printVertexCover();
					auto t_end = std::chrono::high_resolution_clock::now();
					float elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
					runningTimes.push_back(elapsed_time_ms);
				}
				getRunningTime(runningTimes, myfile);
			}
			myfile.close();
		}

		for (int i = 0; i < 4; i++) {
			n = N[i];
			myfile.open("VertexFixed-" + to_string(n) + "-Iterations.csv");
			myfile << "Size(E),";
			myfile << "Mean Time(ms),";
			myfile << "Standard Deviation,";
			myfile << "Standard Error,";
			myfile << "% 90 - CL,";
			myfile << "% 95 - CL\n";

			V = 200;
			E = 200;

			for (int k = 0; k < 10; k++) {
				int new_edge = E + (k * 500);
				myfile << to_string(new_edge) + ",";
				vector<float> runningTimes;

				for (int i = 0; i < n; i++) {
					Graph G(V);
					G.populate(new_edge);
					auto t_start = std::chrono::high_resolution_clock::now();
					G.printVertexCover();
					auto t_end = std::chrono::high_resolution_clock::now();
					float elapsed_time_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();
					runningTimes.push_back(elapsed_time_ms);
				}
				getRunningTime(runningTimes, myfile);
			}
			myfile.close();
		}
	}

	else if (option == 2) {
		int N[] = { 100, 1000, 3000, 5000, 10000 };
		int vertixSize[] = { 5, 6, 7, 8, 9, 10, 15, 20 };


		for (int i = 0; i < 5; i++)
		{
			int n = N[i];
			for (int j = 0; j < 6; j++) {
				double total = 0;
				int vSize = vertixSize[j];
				uniform_int_distribution<int> distribution(0, vSize * (vSize - 1) / 2);
				for (int k = 0; k < n; k++) {
					int randomEdge = distribution(generator);
					Graph g(vSize);
					g.populate(randomEdge);
					//g.printGraph();
					double heuristic = g.printVertexCover();
					//cout << "h: " << heuristic << endl;
					double optimal = g.exactVertexCover();
					//cout << "o: " << optimal << endl;
					total += (optimal / heuristic);
				}
				cout << n << "-" << vSize << ": " << double(total / n) << endl;
			}

		}

	}
	else if (option == 3) {
		vector<int> vertexCover;
		double totalCorrect = 0;
		for (int i = 0; i < 1; i++) {
			uniform_int_distribution<int> distribution1(1000, 10000);
			int vSize = distribution1(generator);
			uniform_int_distribution<int> distribution2(0, 20000);
			int eSize = distribution2(generator);

			Graph G(vSize);
			G.populate(eSize);
			G.returnVertexCover(vertexCover);
			if (G.isVertexCover(vertexCover)) {
				totalCorrect++;
			}

			vertexCover.clear();
		}
		cout << "Correctness Rate: " << totalCorrect / double(1) << endl;
	}

	system("pause");
	return 0;
}