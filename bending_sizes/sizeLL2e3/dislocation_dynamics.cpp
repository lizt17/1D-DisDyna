#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <map>
#include <memory>
#include <iomanip>
#include <cstdlib>
#include <chrono>
#include <sstream>
#include <limits>

#if __cplusplus >= 201703L
  #include <filesystem>
  namespace fs = std::filesystem;
#else
  #include <experimental/filesystem>
  namespace fs = std::experimental::filesystem;
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include "MobilityLaw_W.hpp"
#include "DislocationProperty.hpp"
#include "FileOperation.hpp"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using Clock = std::chrono::high_resolution_clock;

/*=============================================================================*/
/*  dislocation_dynamics_iter3.cpp                                             */
/*-----------------------------------------------------------------------------*/
/*  Main physical idea                                                         */
/*  1) External loading is prescribed by                                       */
/*         Kapp(t) = Kapp0 + KappDot * time                                    */
/*  2) Existing dislocations generate shielding KD and back stress             */
/*  3) Crack-tip effective driving force is                                    */
/*         Ktip = Kapp - KD                                                    */
/*  4) Source resolved shear stress is evaluated by                            */
/*         rss_source = tau_applied + back_stress - tau_image (+ other terms)  */
/*  5) If rss_source >= tau_nuc and Ktip >= Ke, a new dislocation is emitted   */
/*  6) At each nucleation event, write two extra rows:                         */
/*         (a) Ktip just before emission  -> peak point                        */
/*         (b) Ktip just after emission   -> sudden drop point                 */
/*     This makes the Ktip saw-tooth and its peak envelope directly visible    */
/*  7) For each dislocation, compute                                           */
/*         rss_i = tau_applied + tau_interaction - tau_image - tau_friction    */
/*     then get velocity from mobility law                                     */
/*  8) Update positions with stable dt control                                 */
/*  9) Remove dislocations absorbed by crack tip / free surface                */
/* 10) Count how many dislocations have slipped out of the free surface        */
/* 11) Output main history + optional dislocation configuration files          */
/* 12) Stop when Ktip reaches KG_SI or step limit is hit                       */
/*                                                                             */
/*  Pseudocode                                                                 */
/*  -------------------------------------------------------------------------  */
/*  read input.txt                                                             */
/*  set material / geometry / loading / output controls                        */
/*  compute slipLength = LL / cos(theta)                                       */
/*  create output folders                                                      */
/*  initialize dislocation list and event counters                             */
/*  for kInc = 0 ... Nsteps-1                                                  */
/*      sort dislocations by position                                          */
/*      evaluate KD, Ktip, rss_source, rss_i, vel_i                            */
/*      if nucleation criterion is met                                         */
/*          write one event row before nucleation (peak Ktip)                  */
/*          add new dislocation at source                                      */
/*          reevaluate system                                                  */
/*          write one event row after nucleation (drop Ktip)                   */
/*          write one peak row into peak-envelope csv                          */
/*      end if                                                                 */
/*      compute dt from spacing / relative velocity / upper bound              */
/*      move all dislocations                                                  */
/*      count removed-at-tip / removed-at-free-surface                         */
/*      reevaluate end-of-step state                                           */
/*      output regular history row                                             */
/*      if Ktip >= KG_SI: stop                                                 */
/*      time += dt                                                             */
/*  end for                                                                    */
/*=============================================================================*/

namespace term {
static const char* reset  = "\033[0m";
static const char* red    = "\033[31m";
static const char* green  = "\033[32m";
static const char* yellow = "\033[33m";
static const char* blue   = "\033[34m";
static const char* magenta= "\033[35m";
static const char* cyan   = "\033[36m";
}

enum OutputModule {
    MOD_INIT = 0,
    MOD_LOAD = 1,
    MOD_NUC  = 2,
    MOD_MOVE = 3,
    MOD_OUT  = 4,
    MOD_STOP = 5,
    MOD_COUNT = 6
};

enum RecordType {
    REC_REGULAR    = 0,
    REC_NUC_BEFORE = 1,
    REC_NUC_AFTER  = 2,
    REC_FINAL      = 3
};

const char* moduleName(OutputModule mod)
{
    switch (mod)
    {
        case MOD_INIT: return "INIT";
        case MOD_LOAD: return "LOAD";
        case MOD_NUC:  return "NUC";
        case MOD_MOVE: return "MOVE";
        case MOD_OUT:  return "OUT";
        case MOD_STOP: return "STOP";
        default:       return "GEN";
    }
}

const char* moduleColor(OutputModule mod)
{
    switch (mod)
    {
        case MOD_INIT: return term::cyan;
        case MOD_LOAD: return term::blue;
        case MOD_NUC:  return term::green;
        case MOD_MOVE: return term::magenta;
        case MOD_OUT:  return term::yellow;
        case MOD_STOP: return term::red;
        default:       return term::reset;
    }
}

struct InputConfig
{
    /* material */
    double mu_SI;
    double b_SI;
    double rho_SI;
    double nu;
    double tau_nuc_SI;
    double tau_friction_SI;
    double KG_SI;
    double Ke_SI;
    double T;

    /* loading */
    int    loadMode;     // 0 tension, 1 bending
    double KappDot_SI;
    double Kapp0_SI;
    int    useTstress;
    double bending_faw;
    double bending_scale;

    /* geometry */
    double crack_tip;
    double theta;
    double r_source;
    double r_absorbTip;
    double LL;
    double crackLength;

    /* runtime */
    long long Nsteps;
    int outputNum;
    int outputIntervalOverride;
    double shearWaveFraction;
    int maxNumDis;
    double dxMax;

    /* boundary correction */
    int useBoundaryCorrection;
    double bc0;
    double bc1;
    double bc2;
    double bc3;
    double bc4;

    /* output control */
    int writeMainCsv;
    int writeDisConfig;
    int writeSimInfo;
    int writePeakCsv;
    int writeEventRows;
    int writeFinalRow;
    int consoleColor;
    int consoleEnable;
    int module_INIT;
    int module_LOAD;
    int module_NUC;
    int module_MOVE;
    int module_OUT;
    int module_STOP;
    int printEveryStep;
    int printMoveEvery;
    int printOutEvery;
    int flushEveryOut;
};

struct EvalState
{
    double Kapp;
    double KD;
    double Ktip;
    double Klocal;
    double back_stress;
    double rss_source;
    double dxMin;
    double Vmax;

    std::vector<double> currP;
    std::vector<double> tau_int;
    std::vector<double> tau_app;
    std::vector<double> tau_im;
    std::vector<double> rss;
    std::vector<double> currV;
    std::vector<bool>   athermal;
};

struct AdvanceStats
{
    int removeAtTip;
    int removeAtBoundary;
};

class ModuleLogger
{
public:
    ModuleLogger() : useColor_(true), useConsole_(true), writeFile_(true) {}

    void setup(const std::string& simInfoPath, const InputConfig& cfg)
    {
        useColor_   = (cfg.consoleColor != 0);
        useConsole_ = (cfg.consoleEnable != 0);
        writeFile_  = (cfg.writeSimInfo != 0);
        moduleSwitch_[MOD_INIT] = (cfg.module_INIT != 0);
        moduleSwitch_[MOD_LOAD] = (cfg.module_LOAD != 0);
        moduleSwitch_[MOD_NUC ] = (cfg.module_NUC  != 0);
        moduleSwitch_[MOD_MOVE] = (cfg.module_MOVE != 0);
        moduleSwitch_[MOD_OUT ] = (cfg.module_OUT  != 0);
        moduleSwitch_[MOD_STOP] = (cfg.module_STOP != 0);
        if (writeFile_)
        {
            file_.open(simInfoPath.c_str(), std::ios::trunc);
        }
    }

    bool isModuleOn(OutputModule mod) const
    {
        return moduleSwitch_[mod];
    }

    void print(OutputModule mod, const std::string& msg)
    {
        if (!moduleSwitch_[mod]) return;

        std::ostringstream line;
        line << "[" << moduleName(mod) << "] " << msg;

        if (useConsole_)
        {
            if (useColor_)
            {
                std::cout << moduleColor(mod) << line.str() << term::reset << std::endl;
            }
            else
            {
                std::cout << line.str() << std::endl;
            }
        }

        if (writeFile_ && file_.is_open())
        {
            file_ << line.str() << std::endl;
        }
    }

private:
    bool useColor_;
    bool useConsole_;
    bool writeFile_;
    bool moduleSwitch_[MOD_COUNT];
    std::ofstream file_;
};

/*------------------------------------------------------------------------------*/
double getParamOrDefault(const std::map<std::string, double>& params,
                         const std::string& key,
                         double defaultValue)
{
    std::map<std::string, double>::const_iterator it = params.find(key);
    if (it == params.end()) return defaultValue;
    return it->second;
}

long long getLongParamOrDefault(const std::map<std::string, double>& params,
                                const std::string& key,
                                long long defaultValue)
{
    std::map<std::string, double>::const_iterator it = params.find(key);
    if (it == params.end()) return defaultValue;
    return static_cast<long long>(it->second);
}

void ensureDirectory(const std::string& dirName)
{
    fs::create_directories(dirName);
}

InputConfig readInputConfig(const std::string& inputFile)
{
    std::map<std::string, double> params = readFile(inputFile);

    InputConfig cfg;

    cfg.mu_SI             = getParamOrDefault(params, "mu_SI", 161e9);
    cfg.b_SI              = getParamOrDefault(params, "b_SI", 2.74e-10);
    cfg.rho_SI            = getParamOrDefault(params, "rho_SI", 19250.0);
    cfg.nu                = getParamOrDefault(params, "nu", 0.28);
    cfg.tau_nuc_SI        = getParamOrDefault(params, "tau_nuc_SI", 500e6);
    cfg.tau_friction_SI   = getParamOrDefault(params, "tau_friction_SI", 0.0);
    cfg.KG_SI             = getParamOrDefault(params, "KG_SI", 3.1e6);
    cfg.Ke_SI             = getParamOrDefault(params, "Ke_SI", 0.1e6);
    cfg.T                 = getParamOrDefault(params, "T", 400.0);

    cfg.loadMode          = static_cast<int>(getParamOrDefault(params, "loadMode", 0.0));
    cfg.KappDot_SI        = getParamOrDefault(params, "KappDot_SI", 10e6);
    cfg.Kapp0_SI          = getParamOrDefault(params, "Kapp0_SI", 0.0);
    cfg.useTstress        = static_cast<int>(getParamOrDefault(params, "useTstress", 0.0));
    cfg.bending_faw       = getParamOrDefault(params, "bending_faw", 2.0);
    cfg.bending_scale     = getParamOrDefault(params, "bending_scale", 1.0);

    cfg.crack_tip         = getParamOrDefault(params, "crack_tip", 0.0);
    cfg.theta             = getParamOrDefault(params, "theta", M_PI / 4.0);
    cfg.r_source          = getParamOrDefault(params, "r_source", 150.0);
    cfg.r_absorbTip       = getParamOrDefault(params, "r_absorbTip", 5.0);
    cfg.LL                = getParamOrDefault(params, "LL", 1500.0);
    cfg.crackLength       = getParamOrDefault(params, "crackLength", 8e3);

    cfg.Nsteps            = getLongParamOrDefault(params, "Nsteps", static_cast<long long>(5e8));
    cfg.outputNum         = static_cast<int>(getParamOrDefault(params, "outputNum", 500.0));
    cfg.outputIntervalOverride = static_cast<int>(getParamOrDefault(params, "outputIntervalOverride", 0.0));
    cfg.shearWaveFraction = getParamOrDefault(params, "shearWaveFraction", 1e-4);
    cfg.maxNumDis         = static_cast<int>(getParamOrDefault(params, "maxNumDis", 1000.0));
    cfg.dxMax             = getParamOrDefault(params, "dxMax", 100.0);

    cfg.useBoundaryCorrection = static_cast<int>(getParamOrDefault(params, "useBoundaryCorrection", 1.0));
    cfg.bc0               = getParamOrDefault(params, "bc0", 1.0);
    cfg.bc1               = getParamOrDefault(params, "bc1", 0.0);
    cfg.bc2               = getParamOrDefault(params, "bc2", 0.0);
    cfg.bc3               = getParamOrDefault(params, "bc3", 0.0);
    cfg.bc4               = getParamOrDefault(params, "bc4", 0.0);

    cfg.writeMainCsv      = static_cast<int>(getParamOrDefault(params, "writeMainCsv", 1.0));
    cfg.writeDisConfig    = static_cast<int>(getParamOrDefault(params, "writeDisConfig", 1.0));
    cfg.writeSimInfo      = static_cast<int>(getParamOrDefault(params, "writeSimInfo", 1.0));
    cfg.writePeakCsv      = static_cast<int>(getParamOrDefault(params, "writePeakCsv", 1.0));
    cfg.writeEventRows    = static_cast<int>(getParamOrDefault(params, "writeEventRows", 1.0));
    cfg.writeFinalRow     = static_cast<int>(getParamOrDefault(params, "writeFinalRow", 1.0));
    cfg.consoleColor      = static_cast<int>(getParamOrDefault(params, "consoleColor", 1.0));
    cfg.consoleEnable     = static_cast<int>(getParamOrDefault(params, "consoleEnable", 1.0));
    cfg.module_INIT       = static_cast<int>(getParamOrDefault(params, "module_INIT", 1.0));
    cfg.module_LOAD       = static_cast<int>(getParamOrDefault(params, "module_LOAD", 1.0));
    cfg.module_NUC        = static_cast<int>(getParamOrDefault(params, "module_NUC", 1.0));
    cfg.module_MOVE       = static_cast<int>(getParamOrDefault(params, "module_MOVE", 1.0));
    cfg.module_OUT        = static_cast<int>(getParamOrDefault(params, "module_OUT", 1.0));
    cfg.module_STOP       = static_cast<int>(getParamOrDefault(params, "module_STOP", 1.0));
    cfg.printEveryStep    = static_cast<int>(getParamOrDefault(params, "printEveryStep", 0.0));
    cfg.printMoveEvery    = static_cast<int>(getParamOrDefault(params, "printMoveEvery", 0.0));
    cfg.printOutEvery     = static_cast<int>(getParamOrDefault(params, "printOutEvery", 1.0));
    cfg.flushEveryOut     = static_cast<int>(getParamOrDefault(params, "flushEveryOut", 20.0));

    return cfg;
}

/*------------------------------------------------------------------------------*/
inline double cs_SI(const InputConfig& cfg)          { return std::sqrt(cfg.mu_SI / cfg.rho_SI); }
inline double unitSIF(const InputConfig& cfg)        { return cfg.mu_SI * std::sqrt(cfg.b_SI); }
inline double unitSIFrate(const InputConfig& cfg)    { return cfg.mu_SI * cs_SI(cfg) / std::sqrt(cfg.b_SI); }
inline double unitTime(const InputConfig& cfg)       { return cfg.b_SI / cs_SI(cfg); }
inline double mu(const InputConfig&)                 { return 1.0; }
inline double b(const InputConfig&)                  { return 1.0; }
inline double tau_nuc(const InputConfig& cfg)        { return cfg.tau_nuc_SI / cfg.mu_SI; }
inline double tau_friction(const InputConfig& cfg)   { return cfg.tau_friction_SI / cfg.mu_SI; }
inline double Ke(const InputConfig& cfg)             { return cfg.Ke_SI / unitSIF(cfg); }
inline double KappDot(const InputConfig& cfg)        { return cfg.KappDot_SI / unitSIFrate(cfg); }
inline double Kapp0(const InputConfig& cfg)          { return cfg.Kapp0_SI / unitSIF(cfg); }
inline double slipLength(const InputConfig& cfg)
{
    const double c = std::cos(cfg.theta);
    if (std::fabs(c) < 1e-12) return std::numeric_limits<double>::infinity();
    return cfg.LL / std::fabs(c);
}

/*------------------------------------------------------------------------------*/
double schmidFactor(double theta)
{
    return std::sin(theta / 2.0) * std::cos(theta / 2.0) * std::cos(theta / 2.0);
}

double tau_interaction(double ri, double rj, const InputConfig& cfg)
{
    return mu(cfg) * b(cfg) / (2.0 * M_PI * (1.0 - cfg.nu)) / (ri - rj);
}

double tau_image(double r, const InputConfig& cfg)
{
    return mu(cfg) * b(cfg) / (4.0 * M_PI * (1.0 - cfg.nu)) / (r - cfg.crack_tip);
}

double tau_applied(double Kresolved, double r)
{
    if (r <= 0.0) return 0.0;
    return Kresolved / std::sqrt(2.0 * M_PI * r);
}

double tau_Tstress(double Kapp, const InputConfig& cfg)
{
    return Kapp / std::sqrt(M_PI * cfg.crackLength) * 0.4;
}

double tau_bending(double Kapp, double r, const InputConfig& cfg)
{
    double maxS11 = 6.0 * Kapp * 2.0 * std::sqrt(2.0) / std::sqrt(cfg.LL) / cfg.bending_faw;
    return cfg.bending_scale * (maxS11 / std::sqrt(6.0)) * (1.0 - 2.0 * r / slipLength(cfg));
}

double boundary_correction(double r, const InputConfig& cfg)
{
    if (!cfg.useBoundaryCorrection) return 1.0;
    double xi = r / slipLength(cfg);
    if (xi < 0.0) xi = 0.0;
    if (xi > 1.0) xi = 1.0;
    return (((cfg.bc4 * xi + cfg.bc3) * xi + cfg.bc2) * xi + cfg.bc1) * xi + cfg.bc0;
}

double KD_single(double r, const InputConfig& cfg)
{
    if (r <= cfg.crack_tip) return 0.0;
    double pref = 3.0 * (std::sin(cfg.theta) * std::cos(cfg.theta / 2.0));
    return mu(cfg) * b(cfg) / std::sqrt(2.0 * M_PI * r) * pref;
}

double compute_KD(const std::vector<std::shared_ptr<Dislocation> >& disArr,
                  const InputConfig& cfg)
{
    double KD = 0.0;
    for (size_t i = 0; i < disArr.size(); ++i)
    {
        double r = disArr[i]->getPosition();
        KD += boundary_correction(r, cfg) * KD_single(r, cfg);
    }
    return KD;
}

double tau_external_total(double Ktip, double Kapp, double r, const InputConfig& cfg)
{
    // double tau = tau_applied(Ktip * schmidFactor(cfg.theta), r);
    double tau = tau_applied(Kapp * schmidFactor(cfg.theta), r);
    if (cfg.loadMode == 1)
    {
        double tau = tau_applied(Ktip * schmidFactor(cfg.theta), r);
        // tau += tau_bending(Kapp, r, cfg);
    }
    if (cfg.useTstress != 0)
    {
        tau += tau_Tstress(Kapp, cfg);
    }
    return tau;
}

void sortDislocationsDescending(std::vector<std::shared_ptr<Dislocation> >& disArr)
{
    std::sort(disArr.begin(), disArr.end(),
              [](const std::shared_ptr<Dislocation>& a,
                 const std::shared_ptr<Dislocation>& b)
              {
                  return a->getPosition() > b->getPosition();
              });
}

EvalState evaluateSystem(const std::vector<std::shared_ptr<Dislocation> >& disArr,
                         mobilityLaw_W& mobilityLaw,
                         const InputConfig& cfg,
                         double Kapp)
{
    EvalState st;
    int Nd = static_cast<int>(disArr.size());

    st.Kapp        = Kapp;
    st.KD          = compute_KD(disArr, cfg);
    st.Ktip        = st.Kapp - st.KD;
    st.Klocal      = st.Ktip * schmidFactor(cfg.theta);
    st.back_stress = 0.0;
    st.rss_source  = 0.0;
    st.dxMin       = cfg.dxMax;
    st.Vmax        = 1e-12;

    st.currP.resize(Nd);
    st.tau_int.assign(Nd, 0.0);
    st.tau_app.assign(Nd, 0.0);
    st.tau_im.assign(Nd, 0.0);
    st.rss.assign(Nd, 0.0);
    st.currV.assign(Nd, 0.0);
    st.athermal.assign(Nd, false);

    for (int i = 0; i < Nd; ++i)
    {
        st.currP[i] = disArr[i]->getPosition();
        st.back_stress += tau_interaction(cfg.r_source, st.currP[i], cfg);
    }

    st.rss_source = tau_external_total(st.Ktip, st.Kapp, cfg.r_source, cfg)
                  + st.back_stress
                  - tau_image(cfg.r_source, cfg);

    if (Nd == 0)
    {
        st.dxMin = cfg.dxMax;
        return st;
    }

    for (int i = 0; i < Nd - 1; ++i)
    {
        st.dxMin = std::min(st.dxMin, 0.5 * (st.currP[i] - st.currP[i + 1]));
    }
    st.dxMin = std::min(st.dxMin, cfg.dxMax);

    #pragma omp parallel for schedule(static) if(Nd > 1024)
    for (int i = 0; i < Nd; ++i)
    {
        double sum_i = 0.0;
        for (int j = 0; j < Nd; ++j)
        {
            if (j == i) continue;
            sum_i += tau_interaction(st.currP[i], st.currP[j], cfg);
        }
        st.tau_int[i] = sum_i;
    }

    for (int i = 0; i < Nd; ++i)
    {
        // st.tau_app[i] = tau_external_total(st.Ktip, st.Kapp, st.currP[i], cfg);
        st.tau_app[i] = tau_bending(st.Kapp, st.currP[i], cfg); // for bending load
        st.tau_im[i]  = tau_image(st.currP[i], cfg);
        st.rss[i]     = st.tau_app[i] + st.tau_int[i] - st.tau_im[i];

        if (st.rss[i] <= tau_friction(cfg))
        {
            st.rss[i] = 0.0;
        }
        else
        {
            st.rss[i] -= tau_friction(cfg);
        }

        std::pair<bool, double> vel = mobilityLaw.velocity(st.rss[i], cfg.T);
        st.athermal[i] = vel.first;
        st.currV[i] = vel.second;
        st.Vmax = std::max(st.Vmax, std::fabs(st.currV[i]));
    }

    return st;
}

double computePairLimitedDt(const EvalState& st)
{
    int Nd = static_cast<int>(st.currP.size());
    double dtPair = std::numeric_limits<double>::infinity();

    for (int i = 0; i < Nd - 1; ++i)
    {
        double gap = st.currP[i] - st.currP[i + 1];
        double dv  = st.currV[i + 1] - st.currV[i];
        if (dv > 0.0)
        {
            dtPair = std::min(dtPair, 0.45 * gap / dv);
        }
    }
    return dtPair;
}

double computeStableDt(const EvalState& st, const InputConfig& cfg)
{
    double dt = cfg.dxMax / cfg.shearWaveFraction;

    if (!st.currP.empty())
    {
        double vref = std::max(st.Vmax, cfg.shearWaveFraction);
        dt = st.dxMin / vref;

        double dtPair = computePairLimitedDt(st);
        if (std::isfinite(dtPair))
        {
            dt = std::min(dt, dtPair);
        }
    }

    if (!(dt > 0.0) || !std::isfinite(dt))
    {
        dt = cfg.dxMax / cfg.shearWaveFraction;
    }

    return dt;
}

bool shouldNucleate(const EvalState& st, int Nd, const InputConfig& cfg)
{
    return (st.rss_source >= tau_nuc(cfg)) && (st.Ktip >= Ke(cfg)) && (Nd < cfg.maxNumDis);
}

AdvanceStats advanceDislocations(std::vector<std::shared_ptr<Dislocation> >& disArr,
                                 const EvalState& st,
                                 double dt,
                                 const InputConfig& cfg,
                                 int kInc,
                                 ModuleLogger& logger)
{
    std::vector<int> toRemove;
    AdvanceStats adv;
    adv.removeAtTip = 0;
    adv.removeAtBoundary = 0;

    for (size_t i = 0; i < disArr.size(); ++i)
    {
        double newP = st.currP[i] + st.currV[i] * dt;

        disArr[i]->setVelocity(st.currV[i]);
        disArr[i]->setAthermal(st.athermal[i]);
        disArr[i]->setRss(st.rss[i]);

        if (newP <= cfg.r_absorbTip)
        {
            toRemove.push_back(static_cast<int>(i));
            adv.removeAtTip++;
        }
        else if (newP >= slipLength(cfg))
        {
            toRemove.push_back(static_cast<int>(i));
            adv.removeAtBoundary++;
        }
        else
        {
            disArr[i]->setPosition(newP);
        }
    }

    for (int i = static_cast<int>(toRemove.size()) - 1; i >= 0; --i)
    {
        disArr.erase(disArr.begin() + toRemove[i]);
    }

    sortDislocationsDescending(disArr);

    if (logger.isModuleOn(MOD_MOVE) && cfg.printMoveEvery > 0 && (kInc % cfg.printMoveEvery == 0))
    {
        std::ostringstream msg;
        msg << "kInc=" << kInc
            << ", dt=" << dt
            << ", removed_tip=" << adv.removeAtTip
            << ", removed_boundary=" << adv.removeAtBoundary
            << ", Nd_now=" << disArr.size();
        logger.print(MOD_MOVE, msg.str());
    }

    return adv;
}

int getOutputInterval(const InputConfig& cfg)
{
    if (cfg.outputIntervalOverride > 0) return cfg.outputIntervalOverride;
    return static_cast<int>(std::max(1LL, cfg.Nsteps / std::max(1, cfg.outputNum)));
}

void writeMainOutputHeader(std::ofstream& outfile)
{
    outfile << "time,kInc,recordType,Nd,emittedTotal,removedTipStep,removedTipTotal,removedBoundaryStep,removedBoundaryTotal,"
            << "rss_source,Kapp_MPam05,KD_MPam05,Ktip_MPam05,KtipDrop_MPam05,KtipPeak_MPam05,back_stress,dxMin,"
            << "P_first,V_first,RSS_first,P_last,V_last,RSS_last" << std::endl;
}

void writeMainOutputRow(std::ofstream& outfile,
                        double time,
                        long long kInc,
                        int recordType,
                        const EvalState& st,
                        const std::vector<std::shared_ptr<Dislocation> >& disArr,
                        const InputConfig& cfg,
                        int emittedTotal,
                        int removedTipStep,
                        int removedTipTotal,
                        int removedBoundaryStep,
                        int removedBoundaryTotal,
                        double KtipDrop_MPam05,
                        double KtipPeak_MPam05)
{
    outfile << std::scientific << std::setprecision(8)
            << time << ","
            << kInc << ","
            << recordType << ","
            << disArr.size() << ","
            << emittedTotal << ","
            << removedTipStep << ","
            << removedTipTotal << ","
            << removedBoundaryStep << ","
            << removedBoundaryTotal << ","
            << st.rss_source << ","
            << st.Kapp * unitSIF(cfg) / 1e6 << ","
            << st.KD   * unitSIF(cfg) / 1e6 << ","
            << st.Ktip * unitSIF(cfg) / 1e6 << ","
            << KtipDrop_MPam05 << ","
            << KtipPeak_MPam05 << ","
            << st.back_stress << ","
            << st.dxMin;

    if (!disArr.empty())
    {
        outfile << "," << disArr.front()->getPosition()
                << "," << disArr.front()->getVelocity()
                << "," << disArr.front()->getRss()
                << "," << disArr.back()->getPosition()
                << "," << disArr.back()->getVelocity()
                << "," << disArr.back()->getRss();
    }
    else
    {
        outfile << ",0,0,0,0,0,0";
    }
    outfile << std::endl;
}

void writePeakOutputHeader(std::ofstream& peakfile)
{
    peakfile << "time,kInc,emittedTotal,Nd_before,Nd_after,removedBoundaryTotal,"
             << "Kapp_peak_MPam05,Ktip_before_MPam05,Ktip_after_MPam05,Ktip_drop_MPam05"
             << std::endl;
}

void writePeakOutputRow(std::ofstream& peakfile,
                        double time,
                        long long kInc,
                        int emittedTotal,
                        int Nd_before,
                        int Nd_after,
                        int removedBoundaryTotal,
                        double Kapp_peak_MPam05,
                        double Ktip_before_MPam05,
                        double Ktip_after_MPam05)
{
    peakfile << std::scientific << std::setprecision(8)
             << time << ","
             << kInc << ","
             << emittedTotal << ","
             << Nd_before << ","
             << Nd_after << ","
             << removedBoundaryTotal << ","
             << Kapp_peak_MPam05 << ","
             << Ktip_before_MPam05 << ","
             << Ktip_after_MPam05 << ","
             << (Ktip_before_MPam05 - Ktip_after_MPam05) << std::endl;
}

void writeDisConfig(const std::string& disDir,
                    int outFileIndex,
                    const std::vector<std::shared_ptr<Dislocation> >& disArr)
{
    std::ostringstream fileName;
    fileName << disDir << "/evl_" << outFileIndex << ".csv";

    std::ofstream fileDis(fileName.str().c_str(), std::ios::trunc);
    for (size_t i = 0; i < disArr.size(); ++i)
    {
        fileDis << disArr[i]->getId() << ","
                << disArr[i]->isAthermal() << ","
                << std::scientific << std::setprecision(9)
                << disArr[i]->getPosition() << ","
                << disArr[i]->getVelocity() << ","
                << disArr[i]->getRss() << std::endl;
    }
}

/*------------------------------------------------------------------------------*/
int main()
{
    InputConfig cfg = readInputConfig("./input.txt");

    const std::string outputDir = "./output";
    const std::string disConfigDir = "./disConfig";

    ensureDirectory(outputDir);
    ensureDirectory(disConfigDir);

    ModuleLogger logger;
    logger.setup(outputDir + "/simInfo.txt", cfg);

    {
        std::ostringstream msg;
        msg << "Simulation starts. unitSIF=" << unitSIF(cfg)
            << ", unitSIFrate=" << unitSIFrate(cfg)
            << ", unitTime=" << unitTime(cfg);
        logger.print(MOD_INIT, msg.str());
    }
    {
        std::ostringstream msg;
        msg << "loadMode=" << ((cfg.loadMode == 0) ? "tension" : "bending")
            << ", LL=" << cfg.LL
            << ", theta=" << cfg.theta
            << ", slipLength=" << slipLength(cfg)
            << ", r_source=" << cfg.r_source
            << ", useBoundaryCorrection=" << cfg.useBoundaryCorrection;
        logger.print(MOD_INIT, msg.str());
    }

    std::ofstream outfile;
    if (cfg.writeMainCsv != 0)
    {
        std::ostringstream csvName;
        csvName << outputDir << "/outputKD_tau" << static_cast<int>(cfg.tau_nuc_SI / 1e6) << ".csv";
        outfile.open(csvName.str().c_str(), std::ios::trunc);
        writeMainOutputHeader(outfile);
    }

    std::ofstream peakfile;
    if (cfg.writePeakCsv != 0)
    {
        std::ostringstream peakName;
        peakName << outputDir << "/KtipPeaks_tau" << static_cast<int>(cfg.tau_nuc_SI / 1e6) << ".csv";
        peakfile.open(peakName.str().c_str(), std::ios::trunc);
        writePeakOutputHeader(peakfile);
    }

    mobilityLaw_W mobilityLaw;
    std::vector<std::shared_ptr<Dislocation> > disArr;

    int Nd = 0;
    int nextDisID = 0;
    int outFileIndex = 0;
    int removedTipTotal = 0;
    int removedBoundaryTotal = 0;
    double time = 0.0;

    double K0nuc = (tau_nuc(cfg) + tau_friction(cfg) + tau_image(cfg.r_source, cfg))
                 * std::sqrt(2.0 * M_PI * cfg.r_source)
                 / schmidFactor(cfg.theta);
    {
        std::ostringstream msg;
        msg << "Estimated nucleation threshold K0nuc=" << K0nuc * unitSIF(cfg) / 1e6
            << " MPa*m^0.5";
        logger.print(MOD_LOAD, msg.str());
    }

    int outputInterval = getOutputInterval(cfg);

    auto t0 = Clock::now();

    for (long long kInc = 0; kInc < cfg.Nsteps; ++kInc)
    {
        double Kapp = Kapp0(cfg) + KappDot(cfg) * time;

        sortDislocationsDescending(disArr);
        EvalState st = evaluateSystem(disArr, mobilityLaw, cfg, Kapp);
        Nd = static_cast<int>(disArr.size());

        if (logger.isModuleOn(MOD_LOAD) && cfg.printEveryStep > 0 && (kInc % cfg.printEveryStep == 0))
        {
            std::ostringstream msg;
            msg << "kInc=" << kInc
                << ", time=" << time
                << ", Kapp=" << st.Kapp * unitSIF(cfg) / 1e6
                << ", KD=" << st.KD * unitSIF(cfg) / 1e6
                << ", Ktip=" << st.Ktip * unitSIF(cfg) / 1e6
                << ", Nd=" << Nd
                << ", slipOutTotal=" << removedBoundaryTotal;
            logger.print(MOD_LOAD, msg.str());
        }

        if (shouldNucleate(st, Nd, cfg))
        {
            const int Nd_before = Nd;
            const double Kapp_MPam05 = st.Kapp * unitSIF(cfg) / 1e6;
            const double Ktip_before_MPam05 = st.Ktip * unitSIF(cfg) / 1e6;

            if (cfg.writeMainCsv != 0 && cfg.writeEventRows != 0)
            {
                writeMainOutputRow(outfile, time, kInc, REC_NUC_BEFORE, st, disArr, cfg,
                                   nextDisID, 0, removedTipTotal, 0, removedBoundaryTotal,
                                   0.0, Ktip_before_MPam05);
            }

            Dislocation nucleatedDis = {nextDisID, cfg.r_source + b(cfg), 0.0};
            nextDisID++;
            disArr.push_back(std::make_shared<Dislocation>(nucleatedDis));
            sortDislocationsDescending(disArr);
            Nd = static_cast<int>(disArr.size());

            EvalState stAfterNuc = evaluateSystem(disArr, mobilityLaw, cfg, Kapp);
            const double Ktip_after_MPam05 = stAfterNuc.Ktip * unitSIF(cfg) / 1e6;
            const double Ktip_drop_MPam05 = Ktip_before_MPam05 - Ktip_after_MPam05;

            std::ostringstream msg;
            msg << "Nucleating dis id=" << nucleatedDis.getId()
                << ", rss_source=" << st.rss_source
                << ", Kapp=" << Kapp_MPam05
                << ", Ktip(before)=" << Ktip_before_MPam05
                << ", Ktip(after)=" << Ktip_after_MPam05
                << ", drop=" << Ktip_drop_MPam05
                << ", Nd=" << Nd;
            logger.print(MOD_NUC, msg.str());

            if (cfg.writeMainCsv != 0 && cfg.writeEventRows != 0)
            {
                writeMainOutputRow(outfile, time, kInc, REC_NUC_AFTER, stAfterNuc, disArr, cfg,
                                   nextDisID, 0, removedTipTotal, 0, removedBoundaryTotal,
                                   Ktip_drop_MPam05, Ktip_before_MPam05);
            }

            if (cfg.writePeakCsv != 0)
            {
                writePeakOutputRow(peakfile, time, kInc, nextDisID,
                                   Nd_before, Nd, removedBoundaryTotal,
                                   Kapp_MPam05, Ktip_before_MPam05, Ktip_after_MPam05);
            }

            st = stAfterNuc;
        }

        double dt = computeStableDt(st, cfg);
        AdvanceStats adv = advanceDislocations(disArr, st, dt, cfg, static_cast<int>(kInc), logger);
        removedTipTotal += adv.removeAtTip;
        removedBoundaryTotal += adv.removeAtBoundary;

        EvalState stEnd = evaluateSystem(disArr, mobilityLaw, cfg, Kapp);
        Nd = static_cast<int>(disArr.size());

        if (kInc % outputInterval == 0)
        {
            if (cfg.writeMainCsv != 0)
            {
                writeMainOutputRow(outfile, time, kInc, REC_REGULAR, stEnd, disArr, cfg,
                                   nextDisID, adv.removeAtTip, removedTipTotal,
                                   adv.removeAtBoundary, removedBoundaryTotal,
                                   0.0, std::numeric_limits<double>::quiet_NaN());
                if (cfg.flushEveryOut > 0 && (outFileIndex % cfg.flushEveryOut == 0))
                {
                    outfile.flush();
                }
                if (cfg.writePeakCsv != 0 && cfg.flushEveryOut > 0 && (outFileIndex % cfg.flushEveryOut == 0))
                {
                    peakfile.flush();
                }
            }

            if (cfg.writeDisConfig != 0)
            {
                writeDisConfig(disConfigDir, outFileIndex, disArr);
            }

            if (logger.isModuleOn(MOD_OUT) && cfg.printOutEvery > 0 && (outFileIndex % cfg.printOutEvery == 0))
            {
                std::ostringstream msg;
                msg << "outIndex=" << outFileIndex
                    << ", time=" << time
                    << ", Nd=" << Nd
                    << ", Kapp=" << stEnd.Kapp * unitSIF(cfg) / 1e6
                    << ", KD=" << stEnd.KD * unitSIF(cfg) / 1e6
                    << ", Ktip=" << stEnd.Ktip * unitSIF(cfg) / 1e6
                    << ", slipOutStep=" << adv.removeAtBoundary
                    << ", slipOutTotal=" << removedBoundaryTotal;
                logger.print(MOD_OUT, msg.str());
            }

            outFileIndex++;
        }

        if (stEnd.Ktip * unitSIF(cfg) >= cfg.KG_SI)
        {
            std::ostringstream msg;
            msg << "Ktip reached KG. stop at kInc=" << kInc
                << ", time=" << time
                << ", Ktip=" << stEnd.Ktip * unitSIF(cfg) / 1e6
                << ", KG=" << cfg.KG_SI / 1e6
                << ", slipOutTotal=" << removedBoundaryTotal;
            logger.print(MOD_STOP, msg.str());

            if (cfg.writeMainCsv != 0 && cfg.writeFinalRow != 0)
            {
                writeMainOutputRow(outfile, time, kInc, REC_FINAL, stEnd, disArr, cfg,
                                   nextDisID, adv.removeAtTip, removedTipTotal,
                                   adv.removeAtBoundary, removedBoundaryTotal,
                                   0.0, std::numeric_limits<double>::quiet_NaN());
            }
            break;
        }

        time += dt;
    }

    {
        double KappFinal = Kapp0(cfg) + KappDot(cfg) * time;
        EvalState stFinal = evaluateSystem(disArr, mobilityLaw, cfg, KappFinal);
        std::ostringstream msg;
        msg << "Final Nd=" << disArr.size()
            << ", Kapp=" << stFinal.Kapp * unitSIF(cfg) / 1e6
            << ", Ktip=" << stFinal.Ktip * unitSIF(cfg) / 1e6
            << ", emittedTotal=" << nextDisID
            << ", slipOutTotal=" << removedBoundaryTotal;
        logger.print(MOD_STOP, msg.str());
    }
    {
        double elapsed = std::chrono::duration<double>(Clock::now() - t0).count();
        std::ostringstream msg;
        msg << "Wall time=" << elapsed << " s";
        logger.print(MOD_STOP, msg.str());
    }

    if (outfile.is_open()) outfile.close();
    if (peakfile.is_open()) peakfile.close();
    return 0;
}
