#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include "json.hpp"
#include "triangulation.h"

using json = nlohmann::json;
using namespace std;

// Συνάρτηση για φόρτωση των δεδομένων από το JSON αρχείο
void loadDataFromJSON(const string &filename, vector<int> &points_x, vector<int> &points_y, vector<int> &region_boundary, vector<pair<int, int>> &additional_constraints,string &instance_uid)
{
    // Άνοιγμα αρχείου JSON
    ifstream inputFile(filename);
    if (!inputFile.is_open())
    {
        cerr << "Error: Could not open the file " << filename <<    endl;
        return;
    }

    // Φόρτωση JSON δεδομένων
    json j;
    inputFile >> j;

    // Ανάγνωση των arrays από το JSON
    points_x = j["points_x"].get<vector<int>>();
    points_y = j["points_y"].get<vector<int>>();
    region_boundary = j["region_boundary"].get<vector<int>>();
    instance_uid =j["instance_uid"].get<string>();
    cout << "Όνομα που διαβάστηκε: " << instance_uid << endl;
    // Ανάγνωση των constraints
    additional_constraints = j["additional_constraints"].get<vector<pair<int, int>>>();
}
void exportCompletionMessage(string instance_uid) {
    // Δημιουργία ενός JSON αντικειμένου
    cout << instance_uid<<endl;
    json outputData;
    outputData["content_type"] = "CG_SHOP_2025_Solution";
    outputData["instance_uid"] = instance_uid;
    outputData["steiner_points_x"] ="sssss";
    outputData["steiner_points_y"] ="sssss";
    outputData["edges:"] ="sssss";
    

    // Άνοιγμα αρχείου για εγγραφή JSON δεδομένων
    ofstream file("output.json");
    if (file.is_open()) {
        file << outputData.dump(6); // Εγγραφή του JSON με 4 space indentation
        file.close();
        cout << "Το μήνυμα 'complete' αποθηκεύτηκε στο 'output.json'." <<   endl;
        }
         else {
        cerr << "Σφάλμα: Αδυναμία ανοίγματος του αρχείου για εγγραφή." <<   endl;
        }
}


int main()
{
    // Δεδομένα που θα φορτωθούν από το JSON αρχείο
    string instance_uid;
    vector<int> points_x;
    vector<int> points_y;
    vector<int> region_boundary;
    vector<pair<int, int>> additional_constraints;

    // Κάλεσμα της συνάρτησης για φόρτωση δεδομένων
    loadDataFromJSON("data.json", points_x, points_y, region_boundary, additional_constraints,instance_uid);

    // Εκτέλεση τριγωνοποίησης
    triangulate(points_x, points_y, region_boundary, additional_constraints);
    cout << "Όνομα που διαβάστηκε: " << instance_uid << endl;
    exportCompletionMessage(instance_uid);
    return 0;
}
