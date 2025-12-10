#include <iostream>
#include "GeoReadTrack.h"
#include "../BasicEnergySystem/pipe/HPipe.h"

#include <algorithm>

/**
 * GeoReadTrack intended to store an GDAL/OGR geometry and its associated features together
 *
 * Furthermore, it embeds the functionality of the GeoTrack class to work with OGR geometries and features, and provides methods to manipulate the tracks.
 *
 * @brief This class provides methods to access the geometry and the feature. It mainly extends the GeoTrack class to work with OGR geometries and features.
 */
class GeoReadTrack;

// Constructor
GeoReadTrack::GeoReadTrack(OGRGeometry* geom, OGRFeature* feat, const char* diameter_field)
        : geometry(geom), feature(feat) {
    // initialize start and end points from geometry
    init_start_end_points();

    // compute and initialize track length and distance
    GeoReadTrack::compute_track();

    init_diameter(diameter_field);
    init_step_lengths();
}

GeoReadTrack::GeoReadTrack(GeoCoord* start, GeoCoord* end)
        : GeoTrack(start, end), geometry(nullptr), feature(nullptr) {
    // compute and initialize track length and distance
    GeoReadTrack::compute_track();
}

// Copy operator (copies all data but not OGR Pointers)
GeoReadTrack& GeoReadTrack::operator=(const GeoReadTrack& other){
    if (this != &other) {       // not a self-assignment
        start = other.start;
        end = other.end;
        steps = other.steps;
        step_diams = other.step_diams;
        step_lengths = other.step_lengths;
        angles = other.angles;
        diameter = other.diameter;
        length = other.length;
        dist = other.dist;
    }
    return *this;
}

// Destructor
GeoReadTrack::~GeoReadTrack() {
    if (feature)
        OGRFeature::DestroyFeature(feature);
    // geometry pointer is managed by the OGRFeature object and is destroyed when the feature is destroyed
}

// Accessors
OGRGeometry* GeoReadTrack::getGeometry() const { return geometry; }
OGRFeature* GeoReadTrack::getFeature() const { return feature; }

// Methods
void GeoReadTrack::init_diameter(const char* diameter_field_name) {
    diameter = this->getFeature()->GetFieldAsInteger(diameter_field_name);
    step_diams.resize(steps.size()-1);
    fill(step_diams.begin(), step_diams.end(), diameter);
}
void GeoReadTrack::init_step_lengths(double len) {
    if (len < 0)
        len = length;   // default behaviour
    step_lengths.resize(steps.size()-1);
    fill(step_lengths.begin(), step_lengths.end(), len / (double)(steps.size()-1));
}

void GeoReadTrack::get_min_diam(){
    double min_diam = *min_element(step_diams.begin(), step_diams.end());
    if (min_diam <= 0){
        cerr << "Error: Diameter = " << min_diam << endl;
    }
    else
        diameter = min_diam;
}

/**
 * compute_ave_diam() approximates the constant diameter of a pipe that has the same length and the same pressure loss
 * (the Moody factor is assumed to be constant so that L_tot/D_tot^5 = sum(L_i/D_i^5))
 **/
double GeoReadTrack::compute_ave_diam(){
    map<int, double> total_length_of_diam;
    get_total_length();

    for (int i=0; i<step_diams.size(); i++){
        if (step_diams[i] <= 0) {
            cerr << "Error: Diameter = " << step_diams[i] << endl;
            exit(1);
        }
        else
            total_length_of_diam[step_diams[i]] += step_lengths[i];     // a map entry is created and initialized to zero if it does not yet exist
    }

    if (total_length_of_diam.size() == 1)       // all parts are of same diameter
        return step_diams[0];
    else {
        long double denominator = 0.0;
        for (auto &m: total_length_of_diam)
            denominator += m.second * 1e15 / pow(m.first, 5);

        return pow(length / denominator, 0.2);
    }
}

void GeoReadTrack::get_total_length(){
    length = 0.0;
    for (int i=0; i<steps.size()-1; i++){
        if (step_lengths[i] == -1)  // default value
            step_lengths[i] = distance_2d(*steps[i], *steps[i+1]);
        if (step_lengths[i] < 0)
            cerr << "Error Length < 0: Length = " << step_lengths[i] << endl;
        else
            length +=step_lengths[i];
    }
}


/** Method to extract the start and end points coordinates from the geometry as GeoCoordinates
     * The method should handle both LineString and MultiLineString geometries
     * The method assigns the start and end points to the start and end attributes of the GeoTrack object
     * When the geometry is a MultiLineString, the intermediate points are stored in the steps attribute
     */
void GeoReadTrack::init_start_end_points() {
    if (!geometry) {
        throw std::runtime_error("Error: Geometry is null.");
    }
    switch (wkbFlatten(geometry->getGeometryType())) {
        case wkbLineString: {
            OGRLineString* line(geometry->toLineString());
            OGRPoint startPoint, endPoint;
            line->StartPoint(&startPoint);
            line->EndPoint(&endPoint);
            start = new GeoCoord(startPoint.getX(), startPoint.getY()); //TODO coordinate system
            end   = new GeoCoord(endPoint.getX(), endPoint.getY());

            // Keep intermediate points as steps
            steps.clear();
            steps.push_back(start);
            for (int j(1); j < line->getNumPoints() - 1; ++j) {
                OGRPoint point;
                line->getPoint(j, &point);
                steps.push_back(new GeoCoord(point.getX(), point.getY()));
            }
            steps.push_back(end);
            break;
        }
        case wkbMultiLineString: {
            OGRMultiLineString* mline(geometry->toMultiLineString());
            if (mline->getNumGeometries() != 1) {
                std::cerr << "Warning: MultiLineString with more than one geometry, case not fully tested yet." << endl;
            }

            // Keep intermediate points as steps
            for (int i(0); i < mline->getNumGeometries(); ++i) {
                OGRPoint startPoint, endPoint;
                mline->getGeometryRef(i)->StartPoint(&startPoint);
                mline->getGeometryRef(i)->EndPoint(&endPoint);
                start = new GeoCoord(startPoint.getX(), startPoint.getY());
                end = new GeoCoord(endPoint.getX(), endPoint.getY());

                // if (mline->getGeometryRef(i)->getNumPoints() > 2) {
                steps.clear();
                steps.push_back(start);
                for (int j(1); j < mline->getGeometryRef(i)->getNumPoints() - 1; ++j) {
                    OGRPoint point;
                    mline->getGeometryRef(i)->getPoint(j, &point);
                    steps.push_back(new GeoCoord(point.getX(), point.getY()));
                }
                steps.push_back(end);
                // }
            }
            break;
        }
        default: {
            throw std::runtime_error("Error: Geometry type not supported.");
        }
    }
}

/**
 * @brief Overloads the GeoTrack display method to include feature attributes.
 *
 * This method extends the functionality of `GeoTrack::display()` by
 * optionally displaying additional feature attributes when the `isDetailed`
 * parameter is set to `true`.
 *
 * @param isDetailed If true, the method will display detailed feature
 * attributes, otherwise it will display only the basic track information.
 *
 * @see GeoTrack::display()
 */
void GeoReadTrack::display(bool isDetailed) const {
    GeoTrack::display();

    cout << "Diameter: " << diameter << endl;
    cout << "Steps: " << steps.size() << endl;
    cout << "Angles: " << angles.size() << endl;
    cout << endl;

    if(isDetailed){
        std::cout << std::endl;
        std::cout << "Feature ID: " << feature->GetFID() << std::endl;
        std::cout << "Feature attributes: " << std::endl;
        for (int i = 0; i < feature->GetFieldCount(); i++) {
            // print only non-null attributes
            if (feature->IsFieldSet(i)) {
                std::cout << feature->GetFieldDefnRef(i)->GetNameRef() << ": " << feature->GetFieldAsString(i) << std::endl;
            }
        }
        std::cout << "Geometry type: " << OGRGeometryTypeToName(geometry->getGeometryType()) << std::endl;
    }
}

void GeoReadTrack::compute_track(){
    if (start == nullptr || end == nullptr){
        cerr << "Missing definition of start and end of track" << endl;
        exit(1);
    }

    if (dist == -1)
        dist = distance_2d(*start, *end);

    if(steps.size() <= 2) {
        length = dist;
    }
    else { // Compute the length based on the steps
        // Ensure that the start and end points are included in the steps
        if (distance_2d(*start, *steps.front()) > 1e-6) {
            steps.insert(steps.begin(), start);
            // start = steps.front();
        }
        if (distance_2d(*end, *steps.back()) > 1e-6) {
            steps.push_back(end);
            // end = steps.back();
        }
        length = 0;
        for(size_t i(0); i<steps.size()-1; ++i) {
            length += distance_2d(*steps[i], *steps[i + 1]);
        }
    }
}

string GeoReadTrack::write_properties() {
    stringstream ss;
    ss << "\"distance\":" << dist << ",";
    ss << "\"length\":" << length;
    ss << ",\"diameter\":" << diameter;
    // Add angles list if available
    if (!angles.empty()) {
        ss << ",\"angles\":[";
        for (size_t i = 0; i < angles.size(); ++i) {
            ss << angles[i];
            if (i < angles.size() - 1) {
                ss << ",";
            }
        }
        ss << "]";
    }

    string s = std::move(ss).str();
    return s;
}

//Track merging methods

/**
 * Determine the case of merging between two tracks.
 * @param candidate the track to merge with
 * @param tol tolerance for merging
 * @return the case of merging
 */
MergeCase GeoReadTrack::determine_merge_case(const GeoReadTrack* candidate, const double tol) {
    static double tol2, last_tol=0.0;
    double d;

    if (tol != last_tol) {  // recompute tol2 = tol^2 only if tol has changed
        last_tol = tol;
        tol2 = pow(tol, 2);
    }
    if (tol2 > 0) {
        if ((d = this->end->sqr_distance(*candidate->start)) < tol2) { return END_TO_START; }
        if (sqrt(d) - this->length - candidate->length > tol) { return NO_MERGE; }
        if (this->end->sqr_distance(*candidate->end) < tol2) { return END_TO_END; }
        if (this->start->sqr_distance(*candidate->start) < tol2) { return START_TO_START; }
        if (this->start->sqr_distance(*candidate->end) < tol2) { return START_TO_END; }
    }
    else
        cerr << "Merge tolerance is <= 0. Returns without merging." << endl;
    return NO_MERGE;
}


// The merge method could be refactored as a "+" operator overload
// This method could also be adapted to correctly handle attributes and features transmission, for instance:
// - The diameter could be averaged or set to the min value
// - The angles could be concatenated
// These operation might be handled better leveraging the move semantics and the copy constructor of the GeoReadTrack class
/**
 * Merge two tracks based on the determined case
 *
 * Note: The candidate track is left untouched. It must be removed from the vector of tracks if it is to be deleted.
 * @param candidate the track to merge with
 * @param mergeCase the case of merging
 */
void GeoReadTrack::merge_tracks(GeoReadTrack* candidate, const MergeCase mergeCase) {
    switch (mergeCase) {
        case END_TO_START:
            // if order has to be preserved, this end point has to be inserted at the end
            this->steps.insert(this->steps.end(), candidate->steps.begin()+1, candidate->steps.end()); // Append visit's steps: no need to add start point
            this->end = candidate->end; // Update track end

            // append step_diams and step_lengths
            this->step_diams.insert(this->step_diams.end(), candidate->step_diams.begin(), candidate->step_diams.end());
            this->step_lengths.insert(this->step_lengths.end(), candidate->step_lengths.begin(), candidate->step_lengths.end());

            // compute length
            this->length += candidate->length;

            // keep the minimum diameter
            get_min_diam();
            break;

        case END_TO_END:
            // if order has to be preserved, this end point has to be inserted at the end
            this->steps.insert(this->steps.end(), candidate->steps.rbegin()+1, candidate->steps.rend()); // Reverse and append
            this->end = candidate->start; // Update track end

            // append step_diams and step_lengths
            this->step_diams.insert(this->step_diams.end(), candidate->step_diams.rbegin(), candidate->step_diams.rend());
            this->step_lengths.insert(this->step_lengths.end(), candidate->step_lengths.rbegin(), candidate->step_lengths.rend());

            // compute length
            this->length += candidate->length;

            // keep the minimum diameter
            get_min_diam();
            break;

        case START_TO_START:
            // if order has to be preserved, this start point has to be inserted at the beginning
            this->steps.insert(this->steps.begin(), candidate->steps.rbegin(), candidate->steps.rend()-1); // Reverse and prepend
            this->start = candidate->end; // Update track start

            // append step_diams and step_lengths
            this->step_diams.insert(this->step_diams.begin(), candidate->step_diams.rbegin(), candidate->step_diams.rend());
            this->step_lengths.insert(this->step_lengths.begin(), candidate->step_lengths.rbegin(), candidate->step_lengths.rend());

            // compute length
            this->length += candidate->length;

            // keep the minimum diameter
            get_min_diam();
            break;

        case START_TO_END:
            // if order has to be preserved, this start point has to be inserted at the beginning
            this->steps.insert(this->steps.begin(), candidate->steps.begin(), candidate->steps.end()-1); // Prepend steps
            this->start = candidate->start; // Update track start

            // append step_diams and step_lengths
            this->step_diams.insert(this->step_diams.begin(), candidate->step_diams.begin(), candidate->step_diams.end());
            this->step_lengths.insert(this->step_lengths.begin(), candidate->step_lengths.begin(), candidate->step_lengths.end());

            // compute length
            this->length += candidate->length;

            // keep the minimum diameter
            get_min_diam();
            break;

        case NO_MERGE:
            break; // Do nothing
    }
}

/**
 * Detect a tee between this track and a visitor track and returns the iterator of the tee point to be used for splitting
 * @param visitor the track to compare with
 * @param tol tolerance for detecting the tee
 *
 * @return iterator of the tee point
 */
std::vector<GeoCoord *>::iterator GeoReadTrack::detect_tee(const GeoReadTrack *visitor, const double tol) {
    // Iterate through the steps of this track and check visitor's proximity
    for (auto it = steps.begin() + 1; it != steps.end() - 1; ++it) {    // exclude start and end of track
        if (visitor->start->distance(**it) < tol || visitor->end->distance(**it) < tol)
            return it;  // Return the iterator of the first tee found
    }
    return steps.end(); // If no tee is found, return the end iterator
}

bool GeoReadTrack::is_on_track(GeoCoord* coord, const double tol){
    return (position_on_track(coord, tol) != -1);
}

int GeoReadTrack::position_on_track(GeoCoord* coord, const double tol){
    double guess=0.0;
    for (int i=0; i<steps.size(); i++) {
        // if guess < tol, update guess with current distance
        if (guess < tol && (guess = coord->distance(*steps[i])) < tol)
            return i;
        else if (i == steps.size() - 1)     // last step was not a hit
            return -1;
        guess -= step_lengths[i];       // estimate guess
    }
    return -1;
}

set<int> GeoReadTrack::positions_on_track(const vector<GeoCoord*>& coords, double tol, bool check_points){
    set<int> splits;
    GeoCoord* nearest;
    set<GeoCoord*> unused_coords (coords.begin(), coords.end());
    double sqrdist, d, guess=0.0;

    for (int i=0; i<steps.size(); i++) {    // include start and end of track
        // if guess < tol, update guess with current min distance, compare to tol
        if (guess < tol){
            sqrdist = 1e16;
            // update guess with min distance
            for (auto& coord : unused_coords)
                // find minimum sqrdist and corresponding coord
                if ((d=coord->sqr_distance(*steps[i])) < sqrdist) {
                    sqrdist = d;
                    nearest = coord;
                }
            guess = sqrt(sqrdist);

            if (guess < tol) {   // if yes, it is a split, if no, it is a new guess
                splits.insert(i);
                //unused_coords.erase(nearest);     TODO not working
            }
        }
        if (i != steps.size() - 1)
            guess -= step_lengths[i];       // estimate next guess
    }
    if (check_points && (splits.size() != coords.size()))
        throw std::runtime_error("One or more split points not found in track steps!");
    return splits;
}

set<int> GeoReadTrack::positions_on_track(const vector<GeoReadTrack*>& visitors, const double tol){
    set<int> splits;
    //GeoReadTrack* nearest;
    //set<GeoReadTrack*> unused_tracks (visitors.begin(), visitors.end());
    //unused_tracks.erase(this);      // do not compare with itself
    double sqrdist, d, guess=0.0;

    for (int i=0; i<steps.size(); i++) {    // include start and end of track
        // if guess < tol, update guess with current min distance, compare to tol
        if (guess < tol){
            sqrdist = 1e16;
            // update guess with min distance
            for (auto& visitor : visitors) {
                // find minimum sqrdist and corresponding coord
                if (visitor != this && (d = min(visitor->start->sqr_distance(*steps[i]), visitor->end->sqr_distance(*steps[i]))) < sqrdist) {
                    sqrdist = d;
                    //nearest = visitor;
                }
            }
            guess = sqrt(sqrdist);

            if (guess < tol) {   // if yes, it is a split, if no, it is a new guess
                splits.insert(i);
                //unused_tracks.erase(nearest);     //TODO not working
            }
        }
        if (i != steps.size() - 1)
            guess -= step_lengths[i];       // estimate next guess
    }

    return splits;
}

/**
* Method to split the track at multiple given points. For each split point, the track is split into two segments.
* The method takes care of transferring steps from the original track to the new tracks.
* The method assumes the steps are ordered from start to end and base the split on that order.
* Notice that the method doesn't handle the removal of the original track from the vector of tracks and only returns the new tracks in a vector
* @param split_points vector of GeoCoord pointers representing the split points
* @return vector of GeoReadTrack pointers representing the split tracks
*/
std::vector<GeoReadTrack*> GeoReadTrack::split_track(const std::vector<GeoCoord*>& split_points, const double tol, bool check_points) {
    return split_track(positions_on_track(split_points, tol, check_points)); // Return the vector containing all new tracks
}

std::vector<GeoReadTrack*> GeoReadTrack::split_track(const std::vector<GeoReadTrack *> &visitors, double tol){
    return split_track(positions_on_track(visitors, tol));
}

std::vector<GeoReadTrack*> GeoReadTrack::split_track(set<int> split_positions){     // do not call by reference as split_positions is modified
    vector<GeoReadTrack*> new_tracks;
    int old_split = 0;

    // remove 0, add last position
    split_positions.erase(0);
    split_positions.insert((int)steps.size() - 1);

    if (split_positions.size() > 1) {
        for (const auto &split: split_positions) {
            // Create a new track from current_start to split_point
            auto *new_track = new GeoReadTrack(steps[old_split], steps[split]);

            // Copy steps from the original track to the new track up to the split point
            new_track->steps.assign(steps.begin() + old_split, steps.begin() + split + 1);
            new_track->step_diams.assign(step_diams.begin() + old_split, step_diams.begin() + split);
            new_track->step_lengths.assign(step_lengths.begin() + old_split, step_lengths.begin() + split);

            // Get the length for the new track
            new_track->get_total_length();

            // find minimum diameter
            new_track->get_min_diam();

            // Check the start and end points of the new track
            // new_track->check();

            // Add the new track to the split_tracks vector
            new_tracks.push_back(new_track);

            // Update current_start for the next segment
            old_split = split;
        }
    }

    return new_tracks; // Return the vector containing all new tracks
}

/** @overload
  * Method to split the track at a single given point. For each split point, the track is split into two segments.
  * The method takes care of transferring steps from the original track to the new track.
  * Notice that the method doesn't handle the removal of the original track from the vector of tracks and only returns the new tracks in a vector
  * @param split_point single GeoCoord pointer representing the split point
  * @return vector of GeoReadTrack pointers representing the two split tracks
  */
std::vector<GeoReadTrack*> GeoReadTrack::split_track(GeoCoord* split_point, double tol, bool check_points) {
    vector<GeoCoord*> vec_split_points = {split_point};
    return split_track(vec_split_points, tol, check_points);
}

/** Method to extract all tees with a given set of tracks
 * The method relies on the detect_tee method to find the tees between the current track and the tracks in the vector
 * and returns a vector of GeoCoord pointers representing the tees.
 * Note: The detect_tee method assumes that if a tee is found, it is unique and the first one found is returned.
 * @param visitors vector of GeoReadTrack pointers representing the tracks to compare with
 * @return tees vector of GeoCoord pointers representing the tees
 */
std::vector<GeoCoord*> GeoReadTrack::detect_all_tees(const std::vector<GeoReadTrack*>& visitors, double tol) {
    std::vector<GeoCoord*> tees; // Collect all tees for the current track

    for (auto visitor : visitors) {
        if (this != visitor) {
            auto it = this->detect_tee(visitor, tol); // Detect tees between this track and visitor
            if (it != steps.end())
                tees.push_back(*it); // Collect the detected tee point
        }
    }
    return tees; // Return the collected tee points
}

/**
 * Method to calculate the relative angles between the steps of the track
 * The method calculates the relative angles between the steps of the track and stores them in the angles vector.
 * The method uses the cross product to determine the orientation of the angle.
 * The method stores the angles in degrees in the angles vector.
 * @param tol tolerance for storing the angles (in degrees)
 *
 * Note: This method only computes "internal" angles of a track, i.e using angles between three consecutive points in the steps vector.
 * Interfacing angles between different tracks is not handled by this method.
 */
void GeoReadTrack::calculate_relative_angles(double tol) {
    if (steps.size() < 3) {
        return; // No triplets if fewer than 3 points
    }

    for (size_t i = 1; i < steps.size() - 1; ++i) {
        double x1 = steps[i - 1]->east,  y1 = steps[i - 1]->north;
        double x2 = steps[i]->east,      y2 = steps[i]->north;
        double x3 = steps[i + 1]->east,  y3 = steps[i + 1]->north;

        // Vectors
        double v1_x = x2 - x1, v1_y = y2 - y1;
        double v2_x = x3 - x2, v2_y = y3 - y2;

        // Dot product and magnitudes
        double dotProduct = v1_x * v2_x + v1_y * v2_y;
        double mag_v1 = std::sqrt(v1_x * v1_x + v1_y * v1_y);
        double mag_v2 = std::sqrt(v2_x * v2_x + v2_y * v2_y);

        if (mag_v1 == 0 || mag_v2 == 0) {
            continue; // Avoid division by zero
        }

        // Angle (in radians)
        double angle(std::acos(dotProduct / (mag_v1 * mag_v2)));
        // // Cross product to determine orientation
        double crossProduct = v1_x * v2_y - v1_y * v2_x;
        if (crossProduct < 0) {
            angle = -angle; // Make the angle negative for clockwise
        }

        // Convert to degrees
        angle = angle * 180.0 / M_PI;

        // Store angles if bigger than tolerance
        if (std::abs(angle) > tol && std::abs(angle - 180) > tol) {
            angles.push_back(angle);
        }
    }
}

/*** GeoRaccordTrack ***/
void GeoRaccordTrack::get_sst_id(){
    // TODO currently just numbers the tracks by their order
    static int count=1;
    sst_id = count++;
}


/* Outside the class */

/**
 * Loop through a given collection of tracks and merge segments that are close to each other. This handles steps concatenation and track length update.
 * This function only considers start and end points of the tracks and merges them if they are close to each other. Branching/tees is handled separately.
 * Acts directly on the input vector of GeoReadTrack pointers
 * Relies on the determine_merge_case and merge_tracks methods of the GeoReadTrack class.
 * @param tracks vector of GeoReadTrack pointers
 * @param tol tolerance for merging
 */
void merge_all(std::vector<GeoReadTrack*>& tracks, const double tol) {
    bool merged = true;

    while (merged) {
        merged = false; // Reset merge flag

        for (size_t i = 0; i < tracks.size(); ++i) {
            if (tracks[i] != nullptr) {        // Skip already merged tracks
                // Start checking from the next track to avoid self-comparison
                for (size_t j = i + 1; j < tracks.size(); ++j) {
                    if (tracks[j] != nullptr) {        // Skip already merged tracks
                        // Determine which case of merging applies
                        MergeCase mergeCase = tracks[i]->determine_merge_case(tracks[j], tol);

                        if (mergeCase != NO_MERGE) {
                            // Merge tracks based on the determined case
                            tracks[i]->merge_tracks(tracks[j], mergeCase);

                            // remove unused track
                            delete tracks[j];
                            tracks[j] = nullptr;
                            merged = true; // Set merged flag to true
                        }
                    }
                }
            }
        }
        // Remove nullptr from tracks: erase-remove idiom
        tracks.erase(std::remove(tracks.begin(), tracks.end(), nullptr), tracks.end());
    }

    // Recompute the new tracks
    for(auto track : tracks){
        //track->compute_track();
        track->dist = distance_2d(*track->start, *track->end);
        if (track->steps.empty()) {     // should never happen!
            track->get_steps();
            cerr << "Track steps were empty" << endl;
        }
    }
}
/** Method to segment tracks at tees. This allows for the restructuring of the tracks as a tree.
  *
  * Handles the splitting of tracks at tees and the creation of new tracks, as well as the removal of the original tracks when necessary.
  *
  * Directly modify the input vector of GeoReadTrack pointers
  */
void segment_tracks(std::vector<GeoReadTrack*>& tracks, const double tol) {
    segment_tracks(tracks, tracks, tol);
}

void segment_tracks(vector<GeoReadTrack*> &target_tracks, const vector<GeoReadTrack*> &tracks, double tol){
    std::vector<GeoReadTrack*> tracks_to_remove;
    std::vector<GeoReadTrack*> tracks_to_add;

    // Iterate over each track in target_tracks
    for (auto target : target_tracks) {
        std::vector<GeoReadTrack*> new_tracks = target->split_track(tracks, tol);
        if (new_tracks.size() > 1) {    // track was split
            tracks_to_add.insert(tracks_to_add.end(), new_tracks.begin(), new_tracks.end());
            tracks_to_remove.push_back(target);
        }
    }

    // Remove all the original tracks that were split
    for (auto track : tracks_to_remove) {
        // https://en.wikipedia.org/wiki/Erase%E2%80%93remove_idiom
        target_tracks.erase(std::remove(target_tracks.begin(), target_tracks.end(), track), target_tracks.end());
        delete track;
    }

    // Add the newly split tracks to the target_tracks
    target_tracks.insert(target_tracks.end(), tracks_to_add.begin(), tracks_to_add.end());
}