// Datastructures.cc
//
// Student name: Mikko Kaukonen


#include "datastructures.hh"

#include <random>

#include <cmath>

#include <map>
#include <vector>
#include <set>
#include <unordered_set>
#include <iterator>
#include <algorithm>
#include <string>
#include <utility>
#include <iostream>
#include <stack>
#include <list>


// This function was ready in the program code base of the exercise.
// So not done by me!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

std::minstd_rand rand_engine; // Reasonably quick pseudo-random generator

template <typename Type>
Type random_in_range(Type start, Type end)
{
    auto range = end-start;
    ++range;

    auto num = std::uniform_int_distribution<unsigned long int>(0, range-1)(rand_engine);

    return static_cast<Type>(start+num);
}

Datastructures::Datastructures()
{
    // Write any initialization you need here
}

/*!
 * \brief Datastructures::~Datastructures clears the data structures and frees memory.
 * \This destructor releases all allocated memory and clears the data structures.
 */
Datastructures::~Datastructures()
{
    clear_all();
}

/*!
 * \brief Datastructures::Returns the number of stations in the data structure.
 * return The number of stations currently stored in the data structure.
 */
unsigned int Datastructures::station_count()
{
    return stationdata_.size();
}

/*!
 * \brief Datastructures::Clears all data structures and frees allocated memory.
 * This function deallocates memory for all stations and regions stored in the data structures,
 * clears the containers storing them, and clears the trains data structure.
 */
void Datastructures::clear_all()
{
    // Clear stations data structure
    {
        {
        // Iterate through each station and deallocate memory
        for(std::unordered_map<std::string, Station*>::iterator station = stationdata_.begin();
            station != stationdata_.end(); station++)
            {
            delete (station->second);
            }
        // Clear the stations data structure
        stationdata_.clear();
        }
        
        // Clear regions data structure
        {
        // Iterate through each region and deallocate memory
        for(std::unordered_map<unsigned long long int, Region*>::iterator region = regionsdata_.begin();
            region != regionsdata_.end(); region++)
            {
            delete (region->second);
            }
            
        // Clear the regions data structure
        regionsdata_.clear();
        }
        
        // Clear sorting stations data structure
        sortingstation_.clear();

        // Clear trains data structure
        clear_trains();
    }
}

/*!
 * \brief Datastructures::Retrieves all station IDs stored in the data structure.
 * \return A vector containing all station IDs.
 */
std::vector<StationID> Datastructures::all_stations()
{
    // Create a vector to store station names
    std::vector<std::string> station_names;

    // Iterate through each entry in the station data structure
    for(auto & entry :stationdata_)
        // Add the station ID to the vector
        station_names.push_back(entry.first);
    // Return the vector of station IDs
    return station_names;
}

/*!
 * \brief Datastructures::Calculates the distance from the origin (0,0) to a point (x, y)
 * \param x x-coordinate of the point
 * \param y y-coordinate of the point
 * \return The calculated distance as an integer.
 */
int Datastructures::calculate_distance(int x, int y)
{
    // Calculate the distance using the Pythagorean theorem
    int Distance = sqrt(pow(x, 2) + pow(y, 2));
    return Distance;
}

/*!
 * \brief Datastructures::Adds a new station to the data structure 
 * if it doesn't already exist.
 * \param id The ID of the station to be added.
 * \param name The name of the station.
 * \param xy The coordinates of the station.
 * \return rue if the station was successfully added, 
 * false if the station with the same ID already exists.
 */
bool Datastructures::add_station(StationID id, const Name& name, Coord xy)
{
    // Create an empty map for trains associated with the station
    std::map<TrainID, Train*> trains;

    // Check if the station with the given ID already exists
    if(stationdata_.find(id) == stationdata_.end()){

        // If the station doesn't exist, create a new station object
        std::map<unsigned short int,std::vector<std::string>> depatures;

        Station* new_station = new Station{id, name, xy,
                calculate_distance(xy.x, xy.y), xy.y, xy.x, depatures,
                NO_REGION, trains, -1, nullptr, "white", {NO_STATION}, 0 };

        // Insert the new station into the data structures
        stationdata_.insert({id, new_station});
        sortingstation_.insert({name, new_station});
        
        // Return true to indicate that the station was successfully added
        return true;
    }

    else
        // If the station already exists, return false to indicate that it wasn't added
        return false;
}

/*!
 * \brief Datastructures::Retrieves the name of a station given its ID.
 * \param id The ID of the station.
 * \return The name of the station if found, or NO_NAME if the station doesn't exist.
 */
Name Datastructures::get_station_name(StationID id)
{
    // Check if the station with the given ID exists
    if(stationdata_.find(id) != stationdata_.end()){
        // If the station exists, return its name
        return stationdata_[id]->name;
    }
    else
        // If the station doesn't exist, return NO_NAME
        return NO_NAME;
}

/*!
 * \brief Datastructures::etrieves the coordinates of a station given its ID.
 * \param id The ID of the station.
 * \return The coordinates of the station if found, 
 * or NO_COORD if the station doesn't exist.
 */
Coord Datastructures::get_station_coordinates(StationID id)
{
    // Check if the station with the given ID exists
    if(stationdata_.find(id) != stationdata_.end()){
        // If the station exists, return its coordinates
        return stationdata_[id]->coord;
    }
    else
        // If the station doesn't exist, return NO_COORD
        return NO_COORD;
}

/*!
 * \brief Datastructures::Retrieves station IDs sorted alphabetically by station name.
 * \return A vector containing station IDs sorted alphabetically by station name.
 */
std::vector<StationID> Datastructures::stations_alphabetically()
{
    // Create a vector to store sorted station IDs
    std::vector<std::string> sortedstations;
    // Iterate through each entry in the sortingstation_ map
    for(auto & entry : sortingstation_){
        // Add the station ID to the vector
        sortedstations.push_back(entry.second->stationid);
    }

    // Return the vector of sorted station IDs
    return sortedstations;
}

/*!
 * \brief Datastructures::Retrieves station IDs sorted by increasing distance from the origin.
 * \return A vector containing station IDs sorted by increasing distance from the origin.
 */
std::vector<StationID> Datastructures::stations_distance_increasing()
{
    // Create a vector to store pointers to Station objects
    std::vector<Station*> vec;

    // Iterate through each entry in the stationdata_ map
    for(auto & entry : stationdata_){
        // Add the pointer to the Station object to the vector
        vec.push_back(entry.second);
    }
    // Sort the vector of Station pointers by distance from the origin
    std::sort(vec.begin(), vec.end(), [] (Station* a, Station* b)
    {
        // Compare distances between stations
        if(a->Distance != b->Distance)
            {
            return a->Distance < b->Distance;
            }
        // If distances are equal, compare y-coordinates
        else {
            return a->y_coord < b->y_coord;}
            }
    );

    // Create a vector to store sorted station IDs
    std::vector<std::string> sortedcoords;

    // Iterate through the sorted vector of Station pointers
    for(std::vector<Station*>::iterator it = vec.begin(); it != vec.end(); ++it)
        // Add the station ID to the vector
        sortedcoords.push_back((*it)->stationid);

    // Return the vector of sorted station IDs
    return sortedcoords;
}

/*!
 * \brief Datastructures::Finds the station ID associated with the given coordinates.
 * \param xy The coordinates to search for.
 * \return The station ID associated with the given coordinates, or NO_STATION if not found.
 */
StationID Datastructures::find_station_with_coord(Coord xy)
{
    // Iterate through each entry in the stationdata_ map
    for(auto & entry : stationdata_){
        // Check if the station's coordinates match the given coordinates
        if (entry.second->coord == xy)
            // If a station with matching coordinates is found, return its ID
            return entry.second->stationid;
    }
    // If no station with matching coordinates is found, return NO_STATION
    return NO_STATION;
}

/*!
 * \brief Datastructures::Changes the coordinates of a station with the given ID.
 * \param id The ID of the station to update.
 * \param newcoord The new coordinates to assign to the station.
 * \return True if the coordinates were successfully updated, 
 * false if the station with the given ID does not exist.
 */
bool Datastructures::change_station_coord(StationID id, Coord newcoord)
{
    // Check if the station with the given ID exists
    if(stationdata_.find(id) != stationdata_.end()){
        
        // Update the coordinates of the station
        stationdata_[id]->coord=newcoord;

        // Update the distance of the station from the origin
        stationdata_[id]->Distance=(calculate_distance(newcoord.x, newcoord.y));

        // Update the y-coordinate of the station
        stationdata_[id]->y_coord=(newcoord.y);

        // Return true to indicate successful update
        return true;
        }
    else
        // Return false if the station with the given ID does not exist
        return false;
}

/*!
 * \brief Datastructures::Adds a departure time for a train at a specified station.
 * \param stationid The ID of the station.
 * \param trainid The ID of the train.
 * \param time The departure time.
 * \return  True if the departure time for the train was successfully added, 
 * false otherwise.
 */
bool Datastructures::add_departure(StationID stationid, TrainID trainid, Time time)
{
    // Create a vector to store train IDs
    std::vector<TrainID> train;

    // Check if the station with the given ID exists
    if(stationdata_.find(stationid) != stationdata_.end()){

        // Get a reference to the station
        Station* & A = stationdata_.at(stationid);

        // Check if there are no departures at the given time for the station
        if(A->depatures.find(time) == A->depatures.end())
       {
           // Add the train ID to the vector of departures for the given time
            train.push_back(trainid);
            A->depatures.insert({time,train});
           // Return true to indicate successful addition of departure time
            return true;
            }

        // Departures already exist at the given time
        // Check if the train ID already exists in the list of departures for the given time
        else if((std::find(A->depatures[time].begin(),
                            A->depatures[time].end(), trainid) ==
                            A->depatures[time].end())){
                            // If the train ID does not exist, add it to the list of departures for the given time
                            A->depatures[time].push_back(trainid);
            
            // Return true to indicate successful addition of departure time
            return true;
            }
        else
            // Return false to indicate that the departure time for the train already exists
            return false;
            }
    // Station with the given ID does not exist
    else
        // Return false to indicate failure to add departure time
        return false;
}

/*!
 * \brief Datastructures::Removes a departure time for a train at a specified station.
 * \param stationid The ID of the station.
 * \param trainid The ID of the train.
 * \param time The departure time.
 * \return True if the departure time for the train was successfully removed, false otherwise.
 */
bool Datastructures::remove_departure(StationID stationid, TrainID trainid, Time time)
{
    // Check if the station with the given ID exists
    if(stationdata_.find(stationid) != stationdata_.end())
    {
        // Get a reference to the station
        Station* & A = stationdata_[stationid];

        // Find the departure time in the station's departures map
        auto remove_time = A->depatures.find(time);

        // Check if the departure time exists for the station
        if((remove_time != A->depatures.end())
                    // Check if the train ID exists in the list of departures for the given time
                    &(std::find(A->depatures[time].begin(),
                                A->depatures[time].end(), trainid) !=
                                A->depatures[time].end()))
        {
            // Find the iterator pointing to the train ID in the list of departures for the given time
            std::vector<std::string>::iterator remove_train = std::find(A->
                                            depatures[time].begin(), A->
                                            depatures[time].end(), trainid);

            // If the iterator points to the end of the list, the train ID does not exist
            if(remove_train == A->depatures[time].end())
                // Return false to indicate failure to remove departure time
                return false;
            else
                 // Erase the train ID from the list of departures for the given time
                A->depatures[time].erase(remove_train);
            // Return true to indicate successful removal of departure time
            return true;
        }
        // Return true if the departure time does not exist for the station 
        // or if the train ID does not exist in the list of departures for the given time
        return true;

    }
    // Return false if the station with the given ID does not exist
    return false;
}

/*!
 * \brief Comparator function to sort pairs based on departure time and train ID.
 * \param a The first pair to compare.
 * \param b The second pair to compare.
 * \return True if the first pair should precede the second pair in sorting, false otherwise.
 */
bool sortbydepatures(std::pair<Time, TrainID> &a,
              std::pair<Time, TrainID> &b)
{
    // Compare departure times
    if(a.first != b.first)
        // Sort by departure time in ascending order
        return (a.first < b.first);
    else
        // If departure times are the same, sort by train ID in ascending order
        return (a.second < b.second);
}

/*!
 * \brief Datastructures::Retrieves the departures from a station after a specified time.
 * \param stationid The ID of the station.
 * \param time The reference time to filter departures.
 * \return A vector of pairs representing departures after the specified time, sorted by departure time.
 * If the station does not exist or has no departures after the specified time, 
 * a single pair with NO_TIME and NO_TRAIN is returned.
 */
std::vector<std::pair<Time, TrainID>> Datastructures::station_departures_after(StationID stationid, Time time)
{
    // Create a vector to store departures after the specified time
    std::vector<std::pair<Time, TrainID>> vectoreturn;
    vectoreturn.clear();

    // Check if the station with the given ID exists
    if(stationdata_.find(stationid) != stationdata_.end())
    {
        // Get a reference to the station
        Station* & A = stationdata_[stationid];
        // Iterate through the departures of the station
        for(auto & entry : A->depatures)
        {
            // Check if the departure time is after or equal to the specified time
            if(time <= entry.first)
            {
                // Add each departure after the specified time to the return vector
                for(unsigned int i=0; i<entry.second.size(); i++)
                    vectoreturn.push_back(make_pair(entry.first, entry.second.at(i)));
            }
        }
        
        // Sort the vector of departures based on departure time
        std::sort(vectoreturn.begin(), vectoreturn.end(), sortbydepatures);
        // Return the sorted vector of departures after the specified time
        return vectoreturn;
    }
    // If the station does not exist, return a single pair with NO_TIME and NO_TRAIN
    vectoreturn.push_back({NO_TIME, NO_TRAIN});
    return vectoreturn;
}

/*!
 * \brief Adds a new region to the data structure.
 * \param id The ID of the new region.
 * \param name The name of the new region.
 * \param coords The coordinates defining the boundary of the new region.
 * \return True if the region was successfully added, false if a region with the same ID already exists.
 */
bool Datastructures::add_region(RegionID id, const Name &name, std::vector<Coord> coords)
{
    // Check if a region with the given ID already exists
    if(regionsdata_.find(id) == regionsdata_.end())
    {
        // Create vectors to store coordinates, stations, and children of the new region
        std::vector<Coord> coordinates;
        std::vector<StationID> stations;
        std::set<Region*> children;

        // Create and add a new region to the data structure
        Region* new_region = new Region{id, name, coords, stations, children, nullptr};
        regionsdata_.insert({id, new_region});

    // Return true to indicate successful addition of the region
    return true;

    }
    else
        // Return false to indicate a region with the same ID already exists
        return false;
}

/*!
 * \brief Retrieves IDs of all regions in the data structure.
 * \return A vector containing the IDs of all regions.
 */
std::vector<RegionID> Datastructures::all_regions()
{
    // Create a vector to store IDs of all regions
    std::vector<RegionID> regions;

    // Iterate through all entries in the regionsdata_ map
    for(auto & entry: regionsdata_)
        // Retrieve the ID of each region and add it to the vector
        regions.push_back(entry.first);

    // Return the vector containing IDs of all regions
    return regions;
}

/*!
 * \brief Retrieves the name of a region given its ID.
 * \param id The ID of the region.
 * \return The ID of the region.
 */
Name Datastructures::get_region_name(RegionID id)
{
    // Check if the region with the given ID exists
    if(regionsdata_.find(id) != regionsdata_.end())
        // Return the name of the region
        return regionsdata_.at(id)->name;
    else
        // Return NO_NAME if the region does not exist
        return (NO_NAME);
}

/*!
 * \brief Retrieves the coordinates defining a region given its ID.
 * \param id The ID of the region.
 * \return The coordinates defining the boundary of the region if found, 
 * otherwise a vector containing NO_COORD.
 */
std::vector<Coord> Datastructures::get_region_coords(RegionID id)
{
    // Check if the region with the given ID exists
    if(regionsdata_.find(id) != regionsdata_.end())
        // Return the coordinates of the region
        return regionsdata_[id]->coordinates;
    else
        // Return a vector containing NO_COORD if the region does not exist
        return {NO_COORD};
}

/*!
 * \brief Adds a subregion to a parent region in the data structure.
 * \param id The ID of the subregion to be added.
 * \param parentid The ID of the parent region.
 * \return True if the subregion was successfully added to the parent region, false if:
 *         - Either the subregion or the parent region does not exist.
 *         - The subregion already has a parent.
 */
bool Datastructures::add_subregion_to_region(RegionID id, RegionID parentid)
{
    // Get references to the subregion and parent region
    Region* & A = regionsdata_.at(id);
    Region* & B = regionsdata_.at(parentid);

    // Check if the subregion or the parent region does not exist, 
    // or if the subregion already has a parent
    if((regionsdata_.find(id) == regionsdata_.end() ||
        regionsdata_.find(parentid) == regionsdata_.end()) ||
            (regionsdata_[id]->parent != nullptr))
    
        // Return false if any of the conditions are met
        return false;
    else
    {
        // Add the subregion to the parent region
        regionsdata_[id]->parent=B;
        regionsdata_[parentid]->children.insert(A);
        // Return true to indicate successful addition of the subregion to the parent region
        return true;
    }
}

/*!
 * \brief Datastructures::Adds a station to a region in the data structure.
 * \param id The ID of the station to be added.
 * \param parentid The ID of the region to which the station will be added.
 * \return True if the station was successfully added to the region, false if:
 *         - The station or the region does not exist.
 *         - The station already belongs to a region.
 */
bool Datastructures::add_station_to_region(StationID id, RegionID parentid)
{
    // Check if the station or the region does not exist, 
    // or if the station already belongs to a region
    if((stationdata_.find(id) == stationdata_.end()) ||
        (regionsdata_.find(parentid) == regionsdata_.end()) ||
            (stationdata_[id]->region != NO_REGION))
    
        // Return false if any of the conditions are met
        return false;

    else
    {
        // Add the station to the region
        stationdata_[id]->region=parentid;
        regionsdata_[parentid]->stations.push_back(id);

        // Return true to indicate successful addition of the station to the region
        return true;
    }
}

/*!
 * \brief Finds all regions in the hierarchy starting from the given region.
 * \param current Pointer to the current region in the hierarchy.
 * \param vec Vector to store the IDs of the encountered regions.
 * \return Vector containing the IDs of all regions in the hierarchy.
 */
std::vector<RegionID> Datastructures::find_regions(Region* current, std::vector<RegionID> vec)
{
    // Traverse the region hierarchy upwards until reaching the root
    while(current->parent != nullptr)
    {
        // Add the ID of the current region to the vector
        vec.push_back(current->regionid);
        // Move to the parent region
        current = current->parent;
    }
    // Add the ID of the root region to the vector
    vec.push_back(current->regionid);
    // Return the vector containing IDs of all encountered regions
    return vec;
}

/*!
 * \brief Finds all regions containing the given station.
 * \param id The ID of the station.
 * \return Vector containing the IDs of all regions containing the station, 
 * or an empty vector if the station is not associated with any region.
 */
std::vector<RegionID> Datastructures::station_in_regions(StationID id)
{
    // Vector to store the IDs of regions containing the station
    std::vector<RegionID> vec;

    // Check if the station with the given ID exists
    if(stationdata_.find(id) != stationdata_.end())
    {
        // Check if the station is not associated with any region
        if(stationdata_[id]->region == NO_REGION)
            // Return an empty vector if the station is not associated with any region
            return {};
        else
        {
            // Retrieve the current region containing the station
            Region* current = regionsdata_[stationdata_.at(id)->region];

             // Call the find_regions function to find all regions containing the station
            return find_regions(current,vec);
        }
    }
    else
        // Return a vector containing NO_REGION if the station does not exist
        return {NO_REGION};
}

/*!
 * \brief Lists all subregions of the given region.
 * \param region Pointer to the current region whose subregions are to be listed.
 * \param reglist reglist Vector to store the IDs of subregions.
 */
void Datastructures::list_subregions(Region* region, std::vector<RegionID>& reglist)
{
    // Iterator for iterating through subregions
    std::set<Region*> :: iterator reg;
    // Iterate through each subregion of the current region
    for(reg = region->children.begin(); reg != region->children.end(); reg++)
    {
        // Add the ID of the subregion to the vector
        reglist.push_back((*reg)->regionid);
        // Recursively list subregions of the current subregion
        list_subregions(*reg, reglist);
    }

}
/*!
 * \brief Datastructures::Finds all subregions of the region with the specified ID.
 * \param id The ID of the region.
 * \return Vector containing the IDs of all subregions of the specified region, 
 * or an empty vector if the region does not exist or has no subregions.
 */
std::vector<RegionID> Datastructures::all_subregions_of_region(RegionID id)
{
    // Vector to store the IDs of subregions
    std::vector<RegionID> subreg;

    // Check if the region with the given ID exists
    if(regionsdata_.find(id) != regionsdata_.end())
    {
        // Check if the region has no subregions
        if(regionsdata_[id]->children.empty())
            // Return an empty vector if the region has no subregions
            return {};

        else
        {
            // Call the list_subregions function to list all subregions of the specified region
            list_subregions(regionsdata_[id], subreg);
            
            // Return the vector containing IDs of all subregions
            return subreg;
        }


    }
    else
        // Return a vector containing NO_REGION if the region does not exist
        return {NO_REGION};
}

/*!
 * \brief Finds three or less stations closest to the given coordinates.
 * \param xy Coordinates for which the closest stations are to be found.
 * \return Vector containing the IDs of the stations closest to the given coordinates.
 */
std::vector<StationID> Datastructures::stations_closest_to(Coord xy)
{
    // Multimap to store distances and station IDs
    std::multimap<Distance, StationID> closeststations;
    // Vector to store IDs of closest stations
    std::vector<StationID> closestids;

    // Iterate through each station in the data structure
    for(auto & entry: stationdata_)
    {
        // Calculate the distances in x and y dimensions between the station's 
        // coordinates and the given coordinates
        int distx = (entry.second->coord.x - xy.x);
        int disty = (entry.second->y_coord - xy.y);

        // Calculate the distance between the station and the given coordinates
        // and insert the distance and station ID into the multimap
        closeststations.insert({calculate_distance(distx, disty), entry.first});

         // Keep the multimap size limited to 3 (retain only the 3 closest stations)
        if(closeststations.size() > 3)
            closeststations.erase(prev(closeststations.end()));
    }
    
    // Copy the IDs of the closest stations from the multimap to the vector
    for(auto & entry: closeststations)
        closestids.push_back(entry.second);

    // Return the vector containing IDs of the closest stations
    return closestids;
}

/*!
 * \brief Removes the station with the specified ID from the data structure.
 * \param id The ID of the station to be removed.
 * \return True if the station is successfully removed, otherwise false.
 */
bool Datastructures::remove_station(StationID id)
{
    // Check if the station with the given ID exists
    if(stationdata_.find(id) != stationdata_.end())
    {
        // Erase the station's name from the sorting station map
        sortingstation_.erase(stationdata_[id]->name);

        // Find the station in the main station data map
        auto station = stationdata_.find(id);
        
        // Delete the station object
        delete (station->second);

        // Erase the station's entry from the main station data map
        stationdata_.erase(id);

        // Return true to indicate successful removal of the station
        return true;
    }

    else
        // Return false if the station does not exist
        return false;
}

/*!
 * \brief Finds the common parent region of two given regions.
 * \param id1 The ID of the first region.
 * \param id2 The ID of the second region.
 * \return The ID of the common parent region, or NO_REGION 
 * if there is no common parent or one of the regions doesn't exist.
 */
RegionID Datastructures::common_parent_of_regions(RegionID id1, RegionID id2)
{
    // Vector to store the path from region 1 to the root region
    std::vector<RegionID> path1;
    // Vector to store the path from region 2 to the root region
    std::vector<RegionID> path2;
    // Vector to store the common path
    std::vector<RegionID> commonpath;

    // Check if both regions exist in the data structure
    if(regionsdata_.find(id1) != regionsdata_.end() && regionsdata_.find(id2) != regionsdata_.end())
    {
        // Pointer to the first region
        Region* current1 = regionsdata_[id1];
        
        // Find the path from region 1 to the root region
        path1 = find_regions(current1,path1);
        // Pointer to the second region
        Region* current2 = regionsdata_[id2];
        // Find the path from region 2 to the root region
        path2 = find_regions(current2,path2);

        // Iterator for traversing the path of region 1
        auto it1 = path1.size()-1;
        // Iterator for traversing the path of region 2
        auto it2 = path2.size()-1;

        // Compare the paths to find the common parent region
        while(path1[it1] == path2[it2] && it1 >0 && it2 >0)
        {
            commonpath.push_back(path1[it1]);
            it1--;
            it2--;
        }
        // Check if a common parent region exists
        if(commonpath.empty())
            // Return NO_REGION if there is no common parent
            return NO_REGION;

        else
            // Return the ID of the common parent region
            return commonpath.back();
    }
    else
        // Return NO_REGION if one of the specified regions doesn't exist
        return NO_REGION;
}

/*!
 * \brief Datastructures::Adds a train with specified ID and station arrival times.
 * \param trainid The ID of the train to be added.
 * \param stationtimes Vector of pairs containing station IDs and arrival times.
 * \return True if the train is successfully added, otherwise false.
 */
bool Datastructures::add_train(TrainID trainid, std::vector<std::pair<StationID, Time> > stationtimes)
{
    // Vector to store pairs of station IDs and pointers to stations
    std::vector<std::pair<StationID, Station*>> stations;
    // Vector to store station IDs
    std::vector<StationID> list_of_stationids;

    // Check if the train ID already exists in the data structure
    if(traindata_.find(trainid) == traindata_.end())
    {
        // Iterate through each station and its arrival time
        for(auto & entry : stationtimes)
        {
            // Check if the station with the given ID exists
            if(stationdata_.find(entry.first) != stationdata_.end())
            {
                // Add the station ID and pointer to the station to the stations vector
                stations.push_back({entry.first, stationdata_.at(entry.first)});

                // Add the station ID to the list of station IDs
                list_of_stationids.push_back(entry.first);

            }
            else
                // Return false if the station doesn't exist
                return false;
        }
        
        // Create a new Train object with the train ID, station arrival times, stations, and list of station IDs
        Train* new_train = new Train{trainid, stationtimes, stations, list_of_stationids};
        // Insert the new train into the train data map
        traindata_.insert({trainid, new_train});

        // Iterate through each station and its arrival time
        for(auto & entry2 : stationtimes)
        {
            // Update the train data for the station by adding the train ID and pointer to the new train
            stationdata_.at(entry2.first)->trains.insert({trainid, new_train});

            // Add the departure time for the train at the station
            add_departure(entry2.first, trainid, entry2.second);
        }
        // Return true to indicate successful addition of the train
        return true;
    }
    // Return false if the train ID already exists
    return false;
}

/*!
 * \brief Retrieves the IDs of stations reachable from the specified station.
 * \param id The ID of the station to find next reachable stations from.
 * \return A vector containing the IDs of stations reachable from the specified station,
 *         or a vector containing NO_STATION if the specified station does not exist.
 */
std::vector<StationID> Datastructures::next_stations_from(StationID id)
{

    // Vector to store IDs of stations reachable from the specified station
    std::vector<StationID> next_stations;

    // Check if the specified station ID exists in the data structure
    if(stationdata_.find(id) != stationdata_.end())
    {
        // Iterate through each train passing through the specified station
        for(auto & entry : stationdata_[id]->trains)
        {
             // Iterate through each station in the train's route
            for(unsigned long i = 0; i < entry.second->stationids.size() -1; i++)

                // Check if the current station in the route is the specified station
                if(entry.second->stationids[i] == id)

                    // If so, add the ID of the next station in the route to the next_stations vector
                    next_stations.push_back(entry.second->stationids[i+1]);

        }
        
        // Return the vector containing IDs of next reachable stations
        return next_stations;
    }

    // Return a vector containing NO_STATION if the specified station does not exist
    return {NO_STATION};

}

/*!
 * \brief Retrieves the IDs of stations visited by a train after a specified station.
 * \param stationid The ID of the specified station.
 * \param trainid The ID of the specified train.
 * \return A vector containing the IDs of stations visited by the train after the specified station,
 *         or a vector containing NO_STATION if the specified station or train does not exist.
 */
std::vector<StationID> Datastructures::train_stations_from(StationID stationid, TrainID trainid)
{
    // Vector to store IDs of stations visited by the train after the specified station
    std::vector<StationID> stations_after;

    // Check if the specified station ID exists in the data structure 
    // and if the specified train passes through the station
    if(stationdata_.find(stationid) != stationdata_.end() &&
            (stationdata_[stationid]->trains.find(trainid)
             != stationdata_[stationid]->trains.end()))
    {
        // Retrieve the route of the specified train
        stations_after = traindata_[trainid]->stationids;

        // Remove all stations before the specified station from the route
        int i = 0;
        while(stations_after[i] != stationid)
        stations_after.erase(stations_after.begin());

        // Erase the specified station itself from the route
        stations_after.erase(stations_after.begin());

        // If there are stations remaining in the route after removing the specified station, 
        // return them
        if(!stations_after.empty())
            return stations_after;
        else
            
            // If the route is empty after removing the specified station, return NO_STATION
            return {NO_STATION};
}

    // Return NO_STATION if the specified station or train does not exist
    return {NO_STATION};
}

/*!
 * \brief Clears the trains data structure and removes train-related data from stations.
* void Datastructures::clear_trains()
{
    // Clear the trains data structure
    {
    // Iterate through each train in the trains data structure and deallocate memory
    for(std::unordered_map<std::string, Train*>::iterator train = traindata_.begin();
        train != traindata_.end(); train++)
        {
        
        // Deallocate memory for the train
        delete (train->second);
        }
        
    // Remove train-related data from each station
    traindata_.clear();
    }

    for(auto & entry: stationdata_)
        {
        // Clear the departure times for the station
        entry.second->depatures={};
        // Clear the train pointers for the station
        entry.second->trains={};

        }
}

/*!
 * \brief Finds any route between two stations using a depth-first search algorithm.
 * \param fromid The ID of the starting station.
 * \param toid The ID of the destination station.
 * \return A vector of pairs containing station IDs and distances, 
 * representing the route from 'fromid' to 'toid'.
 */
std::vector<std::pair<StationID, Distance>> Datastructures::route_any(StationID fromid, StationID toid)
{
    // Vector to store the route
    std::vector<std::pair<StationID, Distance>> any_route;
    // Stack to store the route from 'fromid' to 'toid'
    std::stack<Station*> route_from_to;

    // Check if both specified stations exist
    if(stationdata_.find(fromid) != stationdata_.end() &&
            stationdata_.find(toid) != stationdata_.end())
    {
        // Set color to white for all stations and initialize other variables
        for(auto it = stationdata_.begin(); it != stationdata_.end(); it++)
        {
            it->second->colour="white";
            it->second->from_node=nullptr;
            it->second->distance_between=0;
            it->second->neighbours={};
        }

       // Set distance of 'fromid' station to 0 and find its neighbors
        stationdata_[fromid]->distance_between = 0;
        stationdata_[fromid]->neighbours = next_stations_from(fromid);

        // Push the 'fromid' station onto the stack and into the return vector
        route_from_to.push(stationdata_[fromid]);
        //any_route.push_back({fromid, stationdata_[fromid]->distance_between});

        // Set the top element of the stack as 'u'
        auto u = route_from_to.top();

        // Repeat until the stack is empty or 'toid' station is reached
        while(!route_from_to.empty() && u!=stationdata_[toid])
        {

            // If 'u' is white, explore its neighbors
            if(u->colour=="white")
            {
                // Change color to grey and explore neighbors
                u->colour="grey";
                auto v = u->neighbours;

                for (auto & element : v)
                {
                    // Calculate distance, set neighbor station and add it to the stack if not already there
                    stationdata_[element]->distance_between = u->distance_between + calculate_distance(
                                    (stationdata_[element]->X_coord -  u->X_coord),
                                    (stationdata_[element]->y_coord -  u->y_coord));
                    stationdata_[element]->from_node=u;
                    stationdata_[element]->neighbours=next_stations_from(element);
                 
                    route_from_to.push(stationdata_[element]);

                    }
                u=route_from_to.top();
                }
            // If 'u' is grey, remove it from the stack
            else if(u->colour=="grey")
            {
                route_from_to.pop();
            }
            // If color is neither white nor grey, return an empty vector
            else
            {
                return {};
            }
        }

        // If destination station is not reached, return an empty vector
        if(u->stationid!=toid)
            return {};
        // If route is found, create the route vector
        else
        {
            // Create a reversed copy of the route vector
            std::vector<std::pair<StationID, Distance>> any_route2;

            while(u!=stationdata_[fromid])
            {
                any_route.push_back({u->stationid, u->distance_between});
                u=u->from_node;
            }
            any_route.push_back({stationdata_[fromid]->stationid,
                                 stationdata_[fromid]->distance_between});

            // Reverse the order of elements in the route vector
            for(auto i=any_route.rbegin(); i<any_route.rend(); i++)
                any_route2.push_back(*i);

            // Return the route vector
            return any_route2;
        }
    }
    // If either of the specified stations doesn't exist, return a vector with NO_STATION
    else
    {
        any_route.push_back({NO_STATION, NO_DISTANCE});
        return any_route;
    }
}

/*!
 * \brief Returns a route with the least number of stations
 * between the specified stations.
 * \param fromid The starting station ID
 * \param toid The destination station ID
 * \return Vector of the IDs and distances of the stations in the route
 */
std::vector<std::pair<StationID, Distance>> Datastructures::route_least_stations(StationID fromid, StationID toid)
{
    // Vector to store the shortest route
    std::vector<std::pair<StationID, Distance>> shortest_route;
    // List to store the shortest route from 'fromid' to 'toid'
    std::list<Station*> shortest_route_from_to;

    // Check if both specified stations exist
    if(stationdata_.find(fromid) != stationdata_.end() &&
            stationdata_.find(toid) != stationdata_.end())
    {
        // Set color to white for all stations and initialize other variables
        for(auto it = stationdata_.begin(); it != stationdata_.end(); it++)
        {
            it->second->colour="white";
            it->second->from_node=nullptr;
            it->second->distance_between=0;
            it->second->neighbours={};
        }

        // Set color to grey for 'fromid' station and find its neighbors
        stationdata_[fromid]->colour="grey";
        stationdata_[fromid]->neighbours = next_stations_from(fromid);
        shortest_route_from_to.push_back(stationdata_[fromid]);

        // Set the front element of the list as 'u'
        auto u = shortest_route_from_to.front();

        // Repeat until the list is empty or 'toid' station is reached
        while (!shortest_route_from_to.empty() && u->stationid!=toid)
        {
            // Remove the front element from the list and change its color to grey
            shortest_route_from_to.pop_front();

            // Explore the neighbors of 'u'
            u->colour="grey";
            auto v = u->neighbours;

            for (auto & element : v)
            {
                // If 'toid' station is reached, calculate the route and return it
                if(stationdata_[element]->stationid==toid)
                {
                    {
                        // Vector to store the route
                        std::vector<std::pair<StationID, Distance>> any_route2;

                        // Calculate distance and set the from node for 'toid' station
                        stationdata_[element]->distance_between = u->distance_between + calculate_distance(
                                        (stationdata_[element]->X_coord -  u->X_coord),
                                        (stationdata_[element]->y_coord -  u->y_coord));
                        stationdata_[element]->from_node=u;

                        // Add 'toid' station to the list
                        shortest_route_from_to.push_back(stationdata_[element]);

                        // Create the route vector
                        while(u!=stationdata_[fromid])
                        {
                            shortest_route.push_back({u->stationid, u->distance_between});
                            u=u->from_node;
                        }
                        shortest_route.push_back({stationdata_[fromid]->stationid,
                                             stationdata_[fromid]->distance_between});

                        // Reverse the order of elements in the route vector
                        for(auto i=shortest_route.rbegin(); i<shortest_route.rend(); i++)
                            any_route2.push_back(*i);

                         // Add 'toid' station and its distance to the route vector
                        any_route2.push_back({toid, stationdata_[toid]->distance_between});

                        // Return the route vector
                        return any_route2;
                    }
                }
                // If neighbor is white, explore it
                else if(stationdata_[element]->colour=="white")
                {
                // Calculate distance, set from node, find neighbors, and change color to grey
                stationdata_[element]->distance_between = u->distance_between + calculate_distance(
                                (stationdata_[element]->X_coord -  u->X_coord),
                                (stationdata_[element]->y_coord -  u->y_coord));
                stationdata_[element]->from_node=u;
                stationdata_[element]->neighbours=next_stations_from(element);
                stationdata_[element]->colour="grey";

                // Add neighbor to the list
                shortest_route_from_to.push_back(stationdata_[element]);

                }
                // If neighbor is grey, change color to black
                else if(stationdata_[element]->colour=="grey")
                    stationdata_[element]->colour="black";


            }
                // Set 'u' as the front element of the list
                u=shortest_route_from_to.front();

            }

        // If destination station is not reached, return an empty vector
        return {};

    }
    
   
    else
    {
         // If either of the specified stations doesn't exist, return a vector with NO_STATION, NO_DISTANCE
        shortest_route.push_back({NO_STATION, NO_DISTANCE});
        return shortest_route;
    }
}

/*!
 * \brief Returns a route that contains a cycle starting from the specified station
 * \param fromid The starting station ID
 * \return Vector of the IDs of the stations in the route
 */
std::vector<StationID> Datastructures::route_with_cycle(StationID fromid)
{
    // Vector to store the IDs of stations in the route
    std::vector<StationID> any_route;
    // Stack to track the traversal of stations in the route
    std::stack<Station*> route_from_to;


    // Check if the station exists
    if(stationdata_.find(fromid) != stationdata_.end() )
    {
        // Reset the stations' attributes
        for(auto it = stationdata_.begin(); it != stationdata_.end(); it++)
        {
            it->second->colour="white";
            it->second->from_node=nullptr;
            it->second->distance_between=0;
            it->second->neighbours={};
        }

        // Get neighbours of the starting station using next_stations_from function
        stationdata_[fromid]->neighbours = next_stations_from(fromid);

        // Add the starting station to the stack
        route_from_to.push(stationdata_[fromid]);

        auto u = route_from_to.top();

        // Loop until the stack is empty
        while(!route_from_to.empty())
        {

            // If the top element of the stack is white
            if(u->colour=="white")
            {
                // Change its color to grey and add its neighbours to the stack
                u->colour="grey";

                auto v = u->neighbours;

                for (auto & element : v)
                {
                    // If the neighbour is white
                    if(stationdata_[element]->colour=="white")
                    {
                        // Store the neighbors and arrival station for each neighbor of u
                        stationdata_[element]->from_node=u;
                        stationdata_[element]->neighbours=next_stations_from(element);

                        // Add the neighbour to the stack
                        route_from_to.push(stationdata_[element]);
                    }
                    else
                    {
                        // If a cycle is detected, construct and return the route
                        std::vector<StationID> any_route2;

                        while(u!=stationdata_[fromid])
                        {
                            // Add current station to the route
                            any_route.push_back(u->stationid);
                            // Move to the previous station
                            u=u->from_node;
                        }
                        
                        // Add the starting station to complete the cycle
                        any_route.push_back(stationdata_[fromid]->stationid);

                        // Reverse the route and add it to any_route2
                        for(auto i=any_route.rbegin(); i<any_route.rend(); i++)
                            any_route2.push_back(*i);
                        
                        // Add the next station after the cycle
                        any_route2.push_back(stationdata_[element]->stationid);

                        // Return the vector containing the route with the cycle
                        return any_route2;

                    }

                    }

                // Remove the top element from the stack, indicating traversal back to a previous station
                u=route_from_to.top();
                }

            
            else
            {
                // Remove the top element from the stack, indicating traversal back to a previous station
                route_from_to.pop();
            }
        }

        // If no cycle is found, return an empty vector
        return {};

    }
   
    else
    {
        // If the starting station does not exist, return a vector containing NO_STATION
        any_route.push_back({NO_STATION});
        return any_route;
    }

}

/*!
 * \brief Finds the shortest distance route between two stations.
 * \param fromid The ID of the starting station.
 * \param toid The ID of the destination station.
 * \return  A vector containing pairs of station IDs 
 * and distances representing the shortest distance route.
 */
std::vector<std::pair<StationID, Distance>> Datastructures::route_shortest_distance(StationID fromid, StationID toid)
{
    // Vector to store the stations and distances of the shortest route.
    std::vector<std::pair<StationID, Distance>> shortest_route;
    // List of stations to traverse to find the shortest route.
    std::list<Station*> shortest_route_from_to;

    // Check if the stations exist
    if(stationdata_.find(fromid) != stationdata_.end() &&
            stationdata_.find(toid) != stationdata_.end())
    {
        // Initialize stations
        for(auto it = stationdata_.begin(); it != stationdata_.end(); it++)
        {
            it->second->colour="white";
            it->second->from_node=nullptr;
            it->second->distance_between=0;
            it->second->neighbours={};
        }

        // Set the color and neighbors of the starting station
        stationdata_[fromid]->colour="grey";
        stationdata_[fromid]->neighbours = next_stations_from(fromid);
        shortest_route_from_to.push_back(stationdata_[fromid]);

        auto u = shortest_route_from_to.front();

        while (!shortest_route_from_to.empty())
        {
            shortest_route_from_to.pop_front();
            u->colour="grey";
            auto v = u->neighbours;

            for (auto & element : v)
            {

                // If the destination station is found
                if(stationdata_[element]->stationid==toid)
                {
                    {
                        // Initialize the vector to store the route
                        std::vector<std::pair<StationID, Distance>> any_route2;

                        // Calculate the distance to the destination station
                        stationdata_[element]->distance_between = u->distance_between + calculate_distance(
                                        (stationdata_[element]->X_coord -  u->X_coord),
                                        (stationdata_[element]->y_coord -  u->y_coord));
                        stationdata_[element]->from_node=u;

                        shortest_route_from_to.push_back(stationdata_[element]);

                        // Store the stations and distances of the route
                        while(u!=stationdata_[fromid])
                        {
                            shortest_route.push_back({u->stationid, u->distance_between});
                            u=u->from_node;
                        }
                        shortest_route.push_back({stationdata_[fromid]->stationid,
                                             stationdata_[fromid]->distance_between});

                        // Reverse the order of the route
                        for(auto i=shortest_route.rbegin(); i<shortest_route.rend(); i++)
                            any_route2.push_back(*i);

                        // Add the destination station and its distance
                        any_route2.push_back({toid, stationdata_[toid]->distance_between});

                        // Return the vector containing the shortest distance route
                        return any_route2;
                    }
                }
                    
                // If the neighboring station is white
                else if(stationdata_[element]->colour=="white")
                {
                // Calculate the distance to the neighboring station
                stationdata_[element]->distance_between = u->distance_between + calculate_distance(
                                (stationdata_[element]->X_coord -  u->X_coord),
                                (stationdata_[element]->y_coord -  u->y_coord));
                stationdata_[element]->from_node=u;
                stationdata_[element]->neighbours=next_stations_from(element);
                stationdata_[element]->colour="grey";

                 // Add the neighboring station to the list of stations to visit
                shortest_route_from_to.push_back(stationdata_[element]);

                }
                // If the neighboring station is grey
                else if(stationdata_[element]->colour=="grey")
                    // If the neighbor has already been visited, mark it as black
                    stationdata_[element]->colour="black";


            }

                u=shortest_route_from_to.front();

            }

        // If neither the destination nor the starting station was found, return an empty vector
        return {};

    }
   
    else
    {
         // If either the starting or destination station was not found, return a vector with NO_STATION and NO_DISTANCE
        shortest_route.push_back({NO_STATION, NO_DISTANCE});
        return shortest_route;
    }
}