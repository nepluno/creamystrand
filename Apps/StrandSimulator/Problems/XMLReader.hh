#ifndef XML_READER
#define XML_READER

#include "../ProblemStepper.hh"
#include "../../../StrandSim/Dynamic/DistanceFieldObject.hh"
#include <rapidxml.hpp>

#include <iostream>
#include <fstream>
#include <limits>
#include <queue>

class XMLReader;

// checkpoint data
struct DumpDataBinary
{
	Scalar m_time;
	int current_frame;
	int current_checkpoint;
	int num_components;
	
	std::string data_fn;
	
	// fluid variables
	VecXx m_x;
	VecXx m_v;
	VecXx m_m;
	VecXx m_radius;
	VecXx m_J;
	VecXx m_vol;
	VecXx m_rest_vol;
	std::vector< int > m_particle_group;
	std::vector< ParticleClassifier > m_classifier;
	VecXuc m_weakened;
	VecXx m_components;
	VecXx m_proj_func;
	
	MatXx m_Fe;
	MatXx m_b;
	MatXx m_B;
	MatXx m_b_trial;
	//	MatXx m_Fe_plus;
	//	MatXx m_b_plus;
	
	// hair variables
	std::vector< VecXx > m_currentDOFs;
	std::vector< VecXx > m_currentAreaDOFs;
	std::vector< VecXx > m_velocities;
	std::vector< VecXx > m_flow_velocities;
	std::vector< VecXx > m_flow_strain;
	std::vector< VecXx > m_flow_components;
	std::vector< Vec2x > m_flow_reservoirs;
	
	// hair goals
	std::vector< std::pair<Vec2i, Vec4x> > m_strand_goals;
	
	// group variables
	std::vector< Vec3x > m_group_pos;
	std::vector< Vec3x > m_group_scale;
	std::vector< Eigen::Quaternion<double> > m_group_rot;
	
	// field variables
	std::vector< int > m_field_i_frames;
	std::vector< Vec3x > m_field_center;
	std::vector< VecXx > m_field_parameter;
	std::vector< Eigen::Quaternion<double> > m_field_rot;
	std::vector< Vec3x > m_field_future_scale;
	
	// mesh variables
	std::vector< int > m_mesh_current_frames;
	std::vector< Vec3xArray > m_previous_mesh_vertices;
	std::vector< Vec3xArray > m_current_mesh_vertices;
	std::vector< Vec3xArray > m_next_mesh_vertices;
	std::vector< Vec3xArray > m_current_mesh_displacement;
};


struct HairPose
{
    double time;
    bool loaded;
    std::vector< VecXx > dofs;
    std::vector< VecXx > areas;
};

struct Script
{
	enum TYPE
	{
		ROTATE,
		TRANSLATE,
        SCALE,
        SWIRL,
        SWITCH,
		
		TYPE_COUNT
	};
	
	enum FUNC
	{
		CUBIC,
		COSINE,
		WENO,
		
		FUNC_COUNT
	};
	
	TYPE type;
	FUNC func;
	int group_index;
	Vec4x v;
	Vec3x origin;
	double start;
	double end;
	double ease_start;
	double ease_end;
	double amplitude;
	double frequency;
	double base_dt;
	double base_pos;
    unsigned int switch_mask;
	
	XMLReader* m_scene;

	bool transform_global;
    bool transform_with_origin;
    bool used;
	
	std::vector<double> base_vertices;
	
	double getNextVelocity( const double& dt, const double& current_time );
	
	void stepScript( const double& dt, const double& current_time );
	
};

class XMLReader : public ProblemStepper
{
public:
    
    XMLReader();
    ~XMLReader();
    
    void PrintAdditionalSettings(const std::string& szfn_dir) override;
	void dumpBinaryCheckpoint(std::string outputdirectory, int current_frame, int current_checkpoint, int file_width) const override;
    int LoadOptions(const char* filename) override;
protected:
	
	template<typename T>
	void loadParam(rapidxml::xml_node<>* nd, const char* name, T& param, std::function<T(const T&)> const& post_process_func = nullptr);

	template<typename T>
	void loadParam(rapidxml::xml_node<>* nd, const char* name, Eigen::Matrix< T, Eigen::Dynamic, 1>& param, int num_components);

	void setupAfterInit(int& current_frame, int& current_check_point) override;
    void setupStrands() override; //TODO: virtual
    void setupMeshes() override; //TODO: virtual
    bool executeScript(int total_substep_id, const Scalar total_substep_dt) override; // execute scripting and any per-step updating for problem
    void setRodCollisionParameters( CollisionParameters& param, const ElasticStrandParameters& strand_param ) override;
    void setSimulationParameters() override;
    
    void projectConstraint(const std::vector<ImplicitStepper*>& steppers) override;

	void loadIntegrator( rapidxml::xml_node<>* node, double& dt );
	void loadParticles(rapidxml::xml_node<>* node, int& maxgroup );
    void loadStrandParameters( rapidxml::xml_node<>* node, const double& dt );
    void loadHairs( rapidxml::xml_node<>* node );
    void loadHairPose( rapidxml::xml_node<>* node );
	void loadCheckpoint( rapidxml::xml_node<>* node, int& current_frame, int& current_check_point );
    void loadCollisionFree( rapidxml::xml_node<>* node );
	void loadScripts( rapidxml::xml_node<>* node);    
    void loadXMLFile( const std::string& filename, std::vector<char>& xmlchars, rapidxml::xml_document<>& doc );
    bool loadTextFileIntoString( const std::string& filename, std::string& filecontents );
    void loadSolidObject( rapidxml::xml_node<>* node, int& maxdf, const double& dt );
    void loadBucketInfo( rapidxml::xml_node<>* node );
    void loadLiquid( rapidxml::xml_node<>* node, const double& dt );
    void initGroups( int num_group );
    void initParticleFaceMapping();
    
    void apply_global_scaling(int group_idx, const Vec3x& mul);
    void apply_local_scaling(int group_idx, const Vec3x& mul);
	void apply_global_rotation(int group_idx, const Eigen::Quaternion<double>& rot);
	void apply_local_rotation(int group_idx, const Eigen::Quaternion<double>& rot);
	void apply_translation(int group_idx, const Vec3x& t);
    void apply_switch(int group_idx, unsigned mask);

	void loadSimpleGravityForces( rapidxml::xml_node<>* node );
    
    void executeHairPose();

    std::string m_filename;
	double m_end_time;
    std::vector<char> m_xmlchars;
    rapidxml::xml_document<> m_doc;
    rapidxml::xml_node<>* m_scene_node;

    std::vector<Vec4x> m_particles;
    std::vector<int> m_particle_groups;
    std::vector<int> m_particle_closed_face_index;
    std::vector<unsigned char> m_fixed;
	std::vector<double> m_surf_height;
	
    std::vector<Scalar> m_solid_flow_height;
    std::vector<DistanceFieldObject> m_fields;
    std::vector<DistanceFieldObject> m_sources;
    std::vector<DistanceFieldObject> m_terminators;
	
    std::vector<Script> m_scripts;
    std::vector<HairPose> m_hairposes;

    std::vector< std::vector< std::pair<int, int> > > m_groups;

    std::vector< std::pair<int, int> > m_global_to_local;
    std::vector< std::vector<int> > m_local_to_global;

    std::vector< std::vector< int > > m_group_fields;
    std::vector< std::vector< int > > m_group_sources;
	std::vector< Vec3x > m_group_pos;
    std::vector< Vec3x > m_group_scale;
	std::vector< Eigen::Quaternion<double> > m_group_rot;
	
	std::vector< Vec3x > m_group_prev_pos;
    std::vector< Vec3x > m_group_prev_scale;
	std::vector< Eigen::Quaternion<double> > m_group_prev_rot;

	Vec3x m_gravity;
    double m_dx;
    double m_bucket_size;
    int m_num_cells;
    
    std::vector< int > m_obj_reader_base;

	friend struct Script;
};

template<typename S, typename T>
inline S lerp_weno( const S value[], T f )
{
	S p1 = value[0] + (value[1] - value[0]) * (f + 2.0) + (value[2] - 2.0 * value[1] + value[0]) * (f + 2.0) * (f + 1.0) * 0.5
	+ (value[3] - 3.0 * value[2] + 3.0 * value[1] - value[0]) * (f + 2.0) * (f + 1.0) * f / 6.0;
	S p2 = value[1] + (value[2] - value[1]) * (f + 1.0) + (value[3] - 2.0 * value[2] + value[1]) * (f + 1.0) * f * 0.5
	+ (value[4] - 3.0 * value[3] + 3.0 * value[2] - value[1]) * (f + 1.0) * f * (f - 1.0) / 6.0;
	S p3 = value[2] + (value[3] - value[2]) * f + (value[4] - 2.0 * value[3] + value[2]) * f * (f - 1.0) * 0.5
	+ (value[5] - 3.0 * value[4] + 3.0 * value[3] - value[2]) * f * (f - 1.0) * (f - 2.0) / 6.0;
	
	T C1 = (2 - f) * (3 - f) / 20.0;
	T C2 = (3 - f) * (f + 2) / 10.0;
	T C3 = (f + 2) * (f + 1) / 20.0;
	
	T IS1 = (814.0 * value[3] * value[3] + 4326 * value[2] * value[2] + 2976 * value[1] * value[1] + 244 * value[0] * value[0] - 3579 * value[2] * value[3] - 6927 * value[2] * value[1]
			 + 1854 * value[2] * value[0] + 2634 * value[3] * value[1] - 683 * value[3] * value[0] - 1659 * value[1] * value[0])   / 180.0;
	T IS2 = (1986 * value[3] * value[3] + 1986 * value[2] * value[2] + 244 * value[1] * value[1] + 244 * value[4] * value[4] + 1074 * value[2] * value[4] - 3777 * value[2] * value[3]
			 - 1269 * value[2] * value[1] + 1074 * value[3] * value[1] - 1269 * value[4] * value[3] - 293 * value[4] * value[1])  / 180.0;
	T IS3 = (814 * value[2] * value[2] + 4326 * value[3] * value[3] + 2976 * value[4] * value[4] + 244 * value[5] * value[5] - 683 * value[2] * value[5] + 2634 * value[2] * value[4]
			 - 3579 * value[2] * value[3] - 6927 * value[3] * value[4] + 1854 * value[3] * value[5] - 1659 * value[4] * value[5]) / 180.0;
	
	const T epsilon = 1e-6;
	T alpha1 = C1 / ((IS1 + epsilon) * (IS1 + epsilon));
	T alpha2 = C2 / ((IS2 + epsilon) * (IS2 + epsilon));
	T alpha3 = C3 / ((IS3 + epsilon) * (IS3 + epsilon));
	
	T sumalpha = alpha1 + alpha2 + alpha3;
	T w1 = alpha1 / sumalpha;
	T w2 = alpha2 / sumalpha;
	T w3 = alpha3 / sumalpha;
	
	return p1 * w1 + p2 * w2 + p3 * w3;
}

inline void max_connected_components(const std::vector<std::vector<int> >& adjs, std::vector< std::vector<int> >& components)
{
	std::vector<bool> visited(adjs.size(), false);
	components.resize(0);

	const int num_verts = (int) adjs.size();
	for(int i = 0; i < num_verts; ++i) {
		if(visited[i]) continue;

		std::vector<int> c;

		std::queue<int> q;
		visited[i] = true;
		q.push(i);

		while(!q.empty()) {
			int u = q.front();
			q.pop();

			c.push_back(u);

			for(int v : adjs[u]) {
				if(visited[v]) continue;
				visited[v] = true;
				q.push(v);
			}
		}

		components.push_back(c);
	}
}

#endif 
