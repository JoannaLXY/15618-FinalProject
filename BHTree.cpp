#include <cstdio>

#include "BHTree.h"

TreeNode::TreeNode(double radius_, double center_x_, double center_y_){
    radius = radius_;
    center.x = center_x_;
    center.y = center_y_;
    // mass = 0.0;
    mass_center = center;
    // num_particles = 0;
}


Quad TreeNode::get_quad(double px, double py){
    if (px <= center.x && py >= center.y){
        return Quad::QuadTopLeft;
    } else if (px >= center.x && py >= center.y){
        return Quad::QuadTopRight;
    } else if (px <= center.x && py <= center.y){
        return Quad::QuadBottomLeft;
    }
    return Quad::QuadBottomRight;
}



void TreeNode::build_subtree(Quad quad){
    double child_x, child_y;
    if (quad == Quad::QuadTopLeft){
        child_x = center.x - radius/2;
        child_y = center.y + radius/2;
    } else if (quad == Quad::QuadTopRight){
        child_x = center.x + radius/2;
        child_y = center.y + radius/2;
    } else if (quad == Quad::QuadBottomLeft){
        child_x = center.x - radius/2;
        child_y = center.y - radius/2;
    } else{
        child_x = center.x + radius/2;
        child_y = center.y - radius/2;
    }
    children[(int) quad] = new TreeNode(radius/2, child_x, child_y);
}


void TreeNode::add_particle(Particle& p){
    if (num_particles == 0){
        num_particles = 1;
        mass = p.mass;
        mass_center = p.position;
        particle = &p;
        return;
    }
    if (num_particles == 1){
        Quad current_quad = get_quad(mass_center.x, mass_center.y);
        build_subtree(current_quad);
        children[(int) current_quad]->add_particle(*particle);
    }
    Quad quad = get_quad(p.position.x, p.position.y);
    if (children[(int) quad] == nullptr){
        build_subtree(quad);
    }
    children[(int) quad]->add_particle(p);

    num_particles++;
    mass_center = (mass_center.mul(mass).add(p.position.mul(p.mass))).mul(1/(mass + p.mass));
    mass += p.mass;
}


void TreeNode::print(){
    printf("mass: %le, mass_x: %le, mass_y: %le, center.x %le, center_y %le, radius %le, num_particles: %d\n", 
            mass, mass_center.x, mass_center.y, center.x, center.y, radius, num_particles);
    for (int i = 0; i < 4; i++){
        if (children[i] != nullptr){
            printf("*******\n");
            printf("Has child %d\n", i);
            children[i]->print();
        }
    }
    printf("-------------\n");
}


vector TreeNode::calculate_force(Particle& p){
    vector r = mass_center.sub(p.position);
    double distance = r.abs();
    if (distance < 1e-20){return r;}    

    if (2 * radius / distance <= THETA || num_particles == 1){
        return r.mul(G * p.mass * mass / pow(distance, 3));
    }
    vector force;
    for (int i = 0; i < 4; i++){
        if (children[i] == nullptr){continue;}
        force = force.add(children[i]->calculate_force(p));
    }
    return force;
}


// TODO: destrcutor
