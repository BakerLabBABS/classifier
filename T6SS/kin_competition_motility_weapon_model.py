
import numpy as np

params = {
    "cell_radius": 0.02,
    "resource_lifetime": 1.0,
    "arena_radius": 1.0,
    "competitor_effective_diffusion": 0.01,
    "background_resource_loss_weight": 0.20,
    "max_relatedness": 0.90,
    "nonmotile_relatedness_length": 1.20,
    "motility_relatedness_decay": 4.00,
    "contact_direct_priority": 0.75,
    "contact_benefit_scale": 0.20,
    "contact_cost": 0.04,
    "diffusible_benefit_scale": 0.70,
    "diffusible_cost": 0.08,
    "competitor_density_default": 1.00,
    "D_nonmotile": 0.01,
    "D_motile": 0.60,
}

def exploitation_length(D, tau=params["resource_lifetime"], a=params["cell_radius"]):
    return np.sqrt(4 * D * tau + a**2)

def relatedness_length(D, xi0=params["nonmotile_relatedness_length"],
                       zeta=params["motility_relatedness_decay"]):
    return xi0 / (1 + zeta * D)

def relatedness_at_distance(r, D, Rmax=params["max_relatedness"]):
    return Rmax * np.exp(-r / relatedness_length(D))

def conspecific_competition_weight(D, lambda_s):
    L = exploitation_length(D)
    return 2 * np.pi * lambda_s * L**2

def focal_access_weight(r, D):
    return np.exp(-r / exploitation_length(D))

def per_kill_inclusive_benefit(r, D, lambda_s):
    w_focal = focal_access_weight(r, D)
    c_conspecific = conspecific_competition_weight(D, lambda_s)
    R = relatedness_at_distance(r, D)
    background = params["background_resource_loss_weight"]
    return (w_focal + R * c_conspecific) / (w_focal + c_conspecific + background)

def contact_per_kill_benefit(D, lambda_s):
    a = params["cell_radius"]
    p = params["contact_direct_priority"]
    return p + (1 - p) * per_kill_inclusive_benefit(a, D, lambda_s)

def contact_encounter_rate(D, lambda_c=params["competitor_density_default"]):
    a = params["cell_radius"]
    R = params["arena_radius"]
    D_rel = D + params["competitor_effective_diffusion"]
    return lambda_c * (4 * np.pi * D_rel) / np.log(R / a)

def diffusible_kill_distance_pdf(r, ell):
    return (r / ell**2) * np.exp(-r / ell)

def diffusible_total_kill_opportunity(ell, lambda_c=params["competitor_density_default"]):
    area = 2 * np.pi * ell**2
    return 1 - np.exp(-lambda_c * area)

def diffusible_mean_per_kill_benefit(D, lambda_s, ell):
    rmax = max(10 * ell, 5 * exploitation_length(D), 3.0)
    r = np.linspace(0, rmax, 2500)
    pdf = diffusible_kill_distance_pdf(r, ell)
    num = np.trapezoid(pdf * per_kill_inclusive_benefit(r, D, lambda_s), r)
    den = np.trapezoid(pdf, r)
    return num / den

def contact_fitness(D, lambda_s, lambda_c=params["competitor_density_default"]):
    return (params["contact_benefit_scale"] *
            contact_encounter_rate(D, lambda_c) *
            contact_per_kill_benefit(D, lambda_s)
            - params["contact_cost"])

def diffusible_fitness(D, lambda_s, ell, lambda_c=params["competitor_density_default"]):
    return (params["diffusible_benefit_scale"] *
            diffusible_total_kill_opportunity(ell, lambda_c) *
            diffusible_mean_per_kill_benefit(D, lambda_s, ell)
            - params["diffusible_cost"])
