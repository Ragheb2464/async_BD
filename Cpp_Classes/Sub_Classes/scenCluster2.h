// Implementation of the KMeans Algorithm
// reference: http://mnemstudio.org/clustering-k-means-example-1.htm



vector <vector <int> > cluster_ids;
vector <vector <int> > cluster_ids_before_equalization;
class Point
{
private:
	int id_point, id_cluster;
	vector<double> values;
	int total_values;
	string name;

public:
	Point(int id_point, vector<double>& values, string name = "")
	{
		this->id_point = id_point;
		total_values = values.size();

		for(int i = 0; i < total_values; i++)
			this->values.push_back(values[i]);

		this->name = name;
		id_cluster = -1;
	}

	int getID()
	{
		return id_point;
	}

	void setCluster(int id_cluster)
	{
		this->id_cluster = id_cluster;
	}

	int getCluster()
	{
		return id_cluster;
	}

	double getValue(int index)
	{
		return values[index];
	}

	int getTotalValues()
	{
		return total_values;
	}

	void addValue(double value)
	{
		values.push_back(value);
	}

	string getName()
	{
		return name;
	}
};

class Cluster
{
private:
	int id_cluster;
	vector<double> central_values;
	vector<Point> points;

public:
	Cluster(int id_cluster, Point point)
	{
		this->id_cluster = id_cluster;

		int total_values = point.getTotalValues();

		for(int i = 0; i < total_values; i++)
			central_values.push_back(point.getValue(i));

		points.push_back(point);
	}

	void addPoint(Point point)
	{
		points.push_back(point);
	}

	bool removePoint(int id_point)
	{
		int total_points = points.size();

		for(int i = 0; i < total_points; i++)
		{
			if(points[i].getID() == id_point)
			{
				points.erase(points.begin() + i);
				return true;
			}
		}
		return false;
	}

	double getCentralValue(int index)
	{
		return central_values[index];
	}

	void setCentralValue(int index, double value)
	{
		central_values[index] = value;
	}

	Point getPoint(int index)
	{
		return points[index];
	}

	int getTotalPoints()
	{
		return points.size();
	}

	int getID()
	{
		return id_cluster;
	}
};

class KMeans
{
private:
	int K; // number of clusters
	int total_values, total_points, max_iterations;
	vector<Cluster> clusters;

	// return ID of nearest center (uses euclidean distance)
	int getIDNearestCenter(Point point)
	{
		double sum = 0.0, min_dist;
		int id_cluster_center = 0;

		for(int i = 0; i < total_values; i++)
		{
			sum += pow(clusters[0].getCentralValue(i) -
					   point.getValue(i), 2.0);
		}

		min_dist = sqrt(sum);

		for(int i = 1; i < K; i++)
		{
			double dist;
			sum = 0.0;

			for(int j = 0; j < total_values; j++)
			{
				sum += pow(clusters[i].getCentralValue(j) -
						   point.getValue(j), 2.0);
			}

			dist = sqrt(sum);

			if(dist < min_dist)
			{
				min_dist = dist;
				id_cluster_center = i;
			}
		}

		return id_cluster_center;
	}

public:
	KMeans(int K, int total_points, int total_values, int max_iterations)
	{
		this->K = K;
		this->total_points = total_points;
		this->total_values = total_values;
		this->max_iterations = max_iterations;
	}

	void run(vector<Point> & points)
	{
		if(K > total_points)
			return;

		vector<int> prohibited_indexes;

		// choose K distinct values for the centers of the clusters
		for(int i = 0; i < K; i++)
		{
			while(true)
			{
				int index_point = rand() % total_points;

				if(find(prohibited_indexes.begin(), prohibited_indexes.end(),
						index_point) == prohibited_indexes.end())
				{
					prohibited_indexes.push_back(index_point);
					points[index_point].setCluster(i);
					Cluster cluster(i, points[index_point]);
					clusters.push_back(cluster);
					break;
				}
			}
		}

		int iter = 1;

		while(true)
		{
			bool done = true;

			// associates each point to the nearest center
			for(int i = 0; i < total_points; i++)
			{
				int id_old_cluster = points[i].getCluster();
				int id_nearest_center = getIDNearestCenter(points[i]);

				if(id_old_cluster != id_nearest_center)
				{
					if(id_old_cluster != -1)
						clusters[id_old_cluster].removePoint(points[i].getID());

					points[i].setCluster(id_nearest_center);
					clusters[id_nearest_center].addPoint(points[i]);
					done = false;
				}
			}
			// recalculating the center of each cluster
			for(int i = 0; i < K; i++)
			{
				for(int j = 0; j < total_values; j++)
				{
					int total_points_cluster = clusters[i].getTotalPoints();
					double sum = 0.0;

					if(total_points_cluster > 0)
					{
						for(int p = 0; p < total_points_cluster; p++)
							sum += clusters[i].getPoint(p).getValue(j);
						clusters[i].setCentralValue(j, sum / total_points_cluster);
					}
				}
			}

			if(done == true || iter >= max_iterations)
			{
				cout << "Break in iteration " << iter << "\n\n";
				break;
			}

			iter++;
		}
		
		
		// write the cluster members before making the equalization moves
		vector <int> aux;
		for(int i = 0; i < K; i++){
			int total_points_cluster =  clusters[i].getTotalPoints();
			aux.clear();
			cout << "Cluster " << clusters[i].getID() + 1 << endl;
			for(int j = 0; j < total_points_cluster; j++){
				aux.push_back(clusters[i].getPoint(j).getID());
				// cout << clusters[i].getPoint(j).getID()  << ": ";
			}
			// cout << endl;
			cluster_ids_before_equalization.push_back(aux);
		}
		
		
		
		
		
		//equalize cardinalities 
		int min_num_scen_cluster = (total_points - (total_points%K))/K;
		while(true)
		{
			//find cluster with max surplus 
			int max_id, max_id_value=0;
			for(int i = 0; i < K; i++)
				if(clusters[i].getTotalPoints() > max_id_value){
					max_id_value = clusters[i].getTotalPoints();
					max_id= i;
				}
			if(max_id_value <= min_num_scen_cluster)
				break;
			//find cluster with min surplus 
			int min_id, min_id_value=total_points+100;
			for(int i = 0; i < K; i++)
				if(clusters[i].getTotalPoints() < min_id_value){
					min_id_value = clusters[i].getTotalPoints();
					min_id= i;
				}
			if(min_id == max_id )
				break;
			bool break_loop=true;
			for(int i=0; i<K; i++)
				if( abs(clusters[i].getTotalPoints() - clusters[max_id].getTotalPoints()) >= 2)
					break_loop=false;
			if(break_loop)
				break;
			//reassign scenarios
			//find the furthest scenario from max_id cluster to be reassigned
			int id_old_cluster = max_id;
			int id_nearest_center = min_id;
			int scen_id = clusters[max_id].getPoint(clusters[max_id].getTotalPoints()-1).getID();//points[clusters[max_id].getTotalPoints()].getID();
			for(int j = 0; j < clusters[max_id].getTotalPoints(); j++)
			{
				double max_dist = 0;
				// int max_dist_id;
				double sum=0;// = pow(clusters[max_id].getCentralValue(j) - points[i].getValue(j), 2.0);
				for(int j = 0; j < total_values; j++)
					sum += pow(clusters[max_id].getCentralValue(j) - points[clusters[max_id].getPoint(j).getID()].getValue(j), 2.0);				
				if(sqrt(sum) > max_dist)
				{
					max_dist=sqrt(sum) ;
					scen_id = clusters[max_id].getPoint(j).getID();
					// cout << endl << scen_id << endl;
				}

			} 
			// cout << endl << "removing scenario " << scen_id << " from cluster " << max_id + 1 << " and adding it to cluster " << min_id +1 << endl;
			// for(int i=0; i < clusters[max_id].getTotalPoints(); i++)
			clusters[id_old_cluster].removePoint(points[scen_id].getID());
			points[scen_id].setCluster(id_nearest_center);
			clusters[id_nearest_center].addPoint(points[scen_id]);	
			for(int i = 0; i < K; i++)//recalaculate the center for each cluster
			{
					for(int j = 0; j < total_values; j++)
					{
						int total_points_cluster = clusters[i].getTotalPoints();
						double sum = 0.0;

						if(total_points_cluster > 0)
						{
							for(int p = 0; p < total_points_cluster; p++)
								sum += clusters[i].getPoint(p).getValue(j);
							clusters[i].setCentralValue(j, sum / total_points_cluster);
						}
					}
			}			
		}
		for(int i = 0; i < K; i++)
			cout << "Cluster-" << i << " ->" << clusters[i].getTotalPoints() <<endl;
		
		// shows elements of clusters
		aux;
		for(int i = 0; i < K; i++)
		{
			int total_points_cluster =  clusters[i].getTotalPoints();

			// cout << "Cluster " << clusters[i].getID() + 1 << endl;
			aux.clear();
			for(int j = 0; j < total_points_cluster; j++)
			{
				// cout << "Point " << clusters[i].getPoint(j).getID()  << ": ";
				aux.push_back(clusters[i].getPoint(j).getID());
				// for(int p = 0; p < total_values; p++)
					// cout << clusters[i].getPoint(j).getValue(p) << " ";

				string point_name = clusters[i].getPoint(j).getName();

				// if(point_name != "")
					// cout << "- " << point_name;

				// cout << endl;
			}			
			cluster_ids.push_back(aux);
			// cout << "Cluster values: ";

			// for(int j = 0; j < total_values; j++)
				// cout << clusters[i].getCentralValue(j) << " ";

			// cout << "\n\n";
		}
	}
};

class ClusMainFunc
{
private:
	int total_points, total_values, K, max_iterations, has_name;
public:
	ClusMainFunc() 
	{ 
	}
	void MainLoop(Data_S *data_S, Search_Param	*search_param, int n_sc, int n_workers)
	{ 
		K = n_workers;//search_param->getNumCluster();
		max_iterations = 1000;
		total_points = n_sc;
		total_values = data_S->getN_od();		
		
		vector<Point> points;
		string point_name;

		for(int i = 0; i < total_points; i++){
			vector<double> values;
			for(int j = 0; j < total_values; j++)
				values.push_back(data_S->getD(j,i));
			point_name = to_string(i);;
			Point p(i, values, point_name);
			points.push_back(p);
		}
		
		KMeans kmeans(K, total_points, total_values, max_iterations);
		
		kmeans.run(points);
		cout <<endl << "I am hereeadadadaddddddddd1111111111111deeeee\n";
	}	
	int GetCluster(int cluster_id, int element) 		{return cluster_ids[cluster_id][element];}
	int GetClusterSize(int cluster_id) 					{return cluster_ids[cluster_id].size();}
	int GetOrginalCluster(int cluster_id, int element) 	{return cluster_ids_before_equalization[cluster_id][element];}
	int GetOriginalClusterSize(int cluster_id) 			{return cluster_ids_before_equalization[cluster_id].size();}
};

