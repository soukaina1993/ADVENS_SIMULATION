// GeoReadTrack.h
#ifndef GEOREADTRACK_H
#define GEOREADTRACK_H

#include <gdal_priv.h>  // GDALDataset
#include "ogrsf_frmts.h"
#include "GeoTrack.h"   // GeoTrack


enum MergeCase {
    NO_MERGE,
    END_TO_START,
    END_TO_END,
    START_TO_START,
    START_TO_END
};

class GeoReadTrack : public GeoTrack {
public:
    // Constructor
    GeoReadTrack()=default;

    GeoReadTrack(OGRGeometry* geom, OGRFeature* feat, const char* diameter_field = "DIAM_NOMIN");

    GeoReadTrack(GeoCoord *start, GeoCoord *end);

    GeoReadTrack(GeoCoord *start, GeoCoord *end, OGRGeometry *geom, OGRFeature *feat);

    // Copy operator (does not copy OGR Pointers)
    GeoReadTrack& operator=(const GeoReadTrack& other);

    // Destructor
    virtual ~GeoReadTrack();

    // Accessors
    OGRGeometry* getGeometry() const;
    OGRFeature* getFeature() const;

    void init_start_end_points();

    void init_diameter(const char* field_name); // Allows to initialize the diameter of the track from a field of the feature
    void init_step_lengths(double len=-1);

    void get_min_diam();        // find minimum diameter in step_diams

    // returns the hydraulically averaged pipe diameter; not yet adapted to available values on pipetable
    double compute_ave_diam();

    void get_total_length();    // sum of all step_lengths

    void display(bool isDetailed = false) const;

    void compute_track() override;

    string write_properties() override;

    // Merging methods
    MergeCase determine_merge_case(const GeoReadTrack* candidate, double tol=0.01);

    void merge_tracks(GeoReadTrack* candidate, MergeCase mergeCase);


    std::vector<GeoCoord *>::iterator detect_tee(const GeoReadTrack *visitor, double tol);
    std::vector<GeoCoord *> detect_all_tees(const std::vector<GeoReadTrack *> &visitors, double tol);

    // Splitting methods
    bool is_on_track(GeoCoord* coord, double tol);
    int position_on_track(GeoCoord* coord, double tol);    // returns idx with steps[idx] == coord, returns -1 if not found
    set<int> positions_on_track(const vector<GeoCoord*> &coords, double tol, bool check_points=false);  // return all split indices
    set<int> positions_on_track(const vector<GeoReadTrack*> &visitors, double tol);

    std::vector<GeoReadTrack*> split_track(set<int> split_positions);   // do not call by reference!
    std::vector<GeoReadTrack*> split_track(GeoCoord* split_point, double tol, bool check_points=false);
    std::vector<GeoReadTrack*> split_track(const std::vector<GeoCoord *> &split_points, double tol, bool check_points=false);
    std::vector<GeoReadTrack*> split_track(const std::vector<GeoReadTrack *> &visitors, double tol);

    void calculate_relative_angles(double tol = 45);

    // Attributes
    vector<double> angles;
    vector<unsigned int> step_diams;    // list of all diameters from start to end
    vector<double> step_lengths;    // list of all lengths from start to end
    unsigned int diameter = 0;

private:
    OGRGeometry* geometry;
    OGRFeature* feature; // IMPORTANT: the logic of the GeoReadTrack points only to the feature of the original OGR object and subsequently added steps are not linked to the feature of theirs (yet?)
};


/*** GeoRaccordTrack ***/
class GeoRaccordTrack : public GeoReadTrack {
public:
    unsigned int sst_id = 0;

    GeoRaccordTrack(OGRGeometry* geom, OGRFeature* feat, const char* diameter_field = "DIAM_NOMIN") : GeoReadTrack(geom, feat, diameter_field) {}
    GeoRaccordTrack(GeoCoord *start, GeoCoord *end) : GeoReadTrack(start, end) {}

    void get_sst_id();      //TODO, but how?

};




/*** Outside the class ***/
void merge_all(std::vector<GeoReadTrack *> &tracks, double tol = 0.01);
// vector<GeoCoord*> extract_all_tees(vector<GeoReadTrack*>& track1, vector<GeoReadTrack*>& track2, bool check_start_end, double tol);
void segment_tracks(vector<GeoReadTrack*> &tracks, double tol);
void segment_tracks(vector<GeoReadTrack*> &target_tracks, const vector<GeoReadTrack*> &tracks, double tol);

inline void check_tracks(std::vector<GeoReadTrack *> &tracks) {
    for (auto track: tracks)
        track->check();
}


#endif // GEOREADTRACK_H

