class Particle{
    public:
        double mass;
        double px, py;
        double vx, vy;
}

enum class Quad{
    QuadTopLeft,
    QuadTopRight,
    QuadBottomLeft,
    QuadBottomRight
}

class TreeNode{
    public:
        double mass;
        double mass_x, mass_y;
        double center_x, center_y;
        double radius;
        TreeNode* children[4];
        int num_particals;

        Quad get_quad(const Particle& p);
        void build_subtree(Quad quad);
        void add_particle(const Particle& p);
        void update_particle(Particle& p);
}