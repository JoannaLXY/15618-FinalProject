#include <cstdio>
#include <unistd.h>
#include <stdlib.h>

#include "BHTree.h"

#define TIME_STEPSIZE 1.0

void write_csv(Particle* p, FILE* fp, int time_idx, int num_particles){
    for (int i = 0; i < num_particles; i++){
        fprintf(fp, "%d,%d,%lf,%lf,%lf,%lf,%lf,%lf\n",
                time_idx, i, time_idx*TIME_STEPSIZE, p[i].position.x, p[i].position.y,
                p[i].velocity.x, p[i].velocity.y, p[i].mass);
    }
}

int main(int argc, char *argv[]){
    int opt = 0;

    char* filename = nullptr;
    char* outfile = nullptr;
    int num_particles;
    int num_iterations;
    double universe_radius = 0.0;

    do {
        opt = getopt(argc, argv, "f:o:i:");
        switch (opt) {
        case 'f':
            filename = optarg;
            break;

        case 'o':
            outfile = optarg;
            break;

        // case 'p':
        //     prob = atof(optarg);
        //     break;

        case 'i':
            num_iterations = atoi(optarg);
            break;

        case -1:
            break;

        default:
            break;
        }
    } while (opt != -1);

    printf("%s\n", filename);

    FILE* fp = fopen(filename,"r");
    fscanf(fp, "%d\n", &num_particles);
    printf("num of particles: %d\n", num_particles);

    fscanf(fp, "%lf\n", & universe_radius);
    printf("universe radius: %le\n", universe_radius);

    Particle* particles = (Particle*) calloc(num_particles, sizeof(Particle));

    char* line;
    size_t len = 0;
    ssize_t read;
    int count = 0;
    while ((read = getline(&line, &len, fp)) != -1) {
        sscanf(line, "%lf %lf %lf %lf %lf %*s", &particles[count].position.x, &particles[count].position.y, 
                &particles[count].velocity.x,  &particles[count].velocity.y, &particles[count].mass);
        count++;
    }
    for (int i= 0; i < num_particles; i++) {
        printf("px: %le, py: %le, vx: %le, vy: %le, mass: %le\n", 
                particles[i].position.x, particles[i].position.y, particles[i].velocity.x, particles[i].velocity.y, particles[i].mass);
    }
    fclose(fp);

    fp = fopen(outfile, "w");
    fprintf(fp, "time_idx,body_idx,t,x,y,vx,vy,m\n");

    write_csv(particles, fp, 0, num_particles);

    for (int iter = 1; iter <= num_iterations; iter++) {
        TreeNode root = TreeNode(universe_radius, 0, 0);
        for (int i = 0; i < num_particles; i++){
            root.add_particle(particles[i]);
        }

        root.print();

        for (int i = 0; i < num_particles; i++){
            vector force = root.calculate_force(particles[i]);
            printf("force ");
            force.print();
            vector acc = force.mul(1.0 / particles[i].mass);
            printf("acc ");
            acc.print();
            vector delta_position = particles[i].velocity.mul(TIME_STEPSIZE)
                .add(acc.mul(0.5*TIME_STEPSIZE*TIME_STEPSIZE));
            printf("delta_x ");
            delta_position.print();
            particles[i].position = particles[i].position.add(delta_position);
            particles[i].velocity = particles[i].velocity.add(acc.mul(TIME_STEPSIZE));
        }
        write_csv(particles, fp, iter, num_particles);
    }

    fclose(fp);

    free(particles);

    return 0;
}