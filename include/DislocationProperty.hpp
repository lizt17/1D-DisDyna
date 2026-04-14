#ifndef DISLOCATIONPROPERTY_H
#define DISLOCATIONPROPERTY_H

class Dislocation {
public:
    // Constructor
    Dislocation(int id, double position, double velocity, bool athermal = false)
        : id(id), position(position), velocity(velocity) {}
    // Getters
    int getId() const { return id; }
    double getPosition() const { return position; }
    double getVelocity() const { return velocity; }
    bool isAthermal() const { return athermal; }
    double getRss() const { return rss; }

    // Setters
    void setPosition(double newPosition) 
    {
        position = newPosition;
    }

    void setVelocity(double newVelocity) 
    {
        velocity = newVelocity;
    }

    void setAthermal(bool isAthermal) 
    {
        athermal = isAthermal;
    }

    void setRss(double newRss) 
    {
        rss = newRss;
    }

private:
    int id;          // Unique ID for the dislocation
    double position; // Position in 1D
    double velocity; // Velocity in 1D
    bool athermal; // Athermal state
    double rss; // Resolved shear stress
};


#endif // DISLOCATIONPROPERTY_H
