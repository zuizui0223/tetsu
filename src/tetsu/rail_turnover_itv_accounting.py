"""Exact partition of community-weighted mean (CWM) change into species
turnover versus within-species trait deviation with an external taxon reference.

Pure mathematical preparatory tool; deliberately does not derive community
abundance, flower dates, height, optimal traits or railway preference from
the iNaturalist photo opportunity counts. Real original rail-vs-matched
nonrail fixed-area community weights, within-species trait means in BOTH
habitats and independent reference values must be provided first.

For species i let p_iR and p_iC be normalized nonnegative observed weights
and z_iR and z_iC true observed within-site mean traits. For external
z_i* shared in both cases:
  ΔCWM = Σ (p_iR - p_iC) z_i*
       + Σ p_iR(z_iR - z_i*) - Σ p_iC(z_iC - z_i*)

Composition and ITV terms depend on chosen external reference and source
sampling design. This is an EXACT accounting identity, NOT causation,
adaptive fitness or a demonstration of heritable local adaptation.
Missing reference value for ANY observed species is not quietly imputed,
and no photo-based phenophase absent annotation is scored as zero flower.
"""
from __future__ import annotations
import math
from collections.abc import Mapping

def _prepare(habitat,kind):
    if not isinstance(habitat,Mapping) or not habitat:
        raise ValueError(kind+" must be an observed nonempty species->(weight,trait) mapping")
    total=0.
    rows={}
    for species,record in habitat.items():
        if not isinstance(species,str) or not species.strip():
            raise ValueError("All observed species must have nonblank exact source taxon names")
        if not isinstance(record,(tuple,list)) or len(record)!=2:
            raise ValueError(kind+" observed species need numeric weight and trait mean")
        w,z=record
        if not all(isinstance(x,(int,float)) and not isinstance(x,bool) and
                   math.isfinite(x) for x in (w,z)) or w<=0:
            raise ValueError(kind+" must use observed positive finite weights and trait means; absence is omitted")
        rows[species]=(float(w),float(z))
        total+=w
    if total<=0 or not math.isfinite(total):
        raise ValueError(kind+" must have finite positive total weight")
    return {tax:(w/total,z) for tax,(w,z) in rows.items()}

def decompose_rail_vs_control(rail,control,external_reference):
    r=_prepare(rail,"rail")
    c=_prepare(control,"control")
    if not isinstance(external_reference,Mapping):
        raise ValueError("Independent species-level reference traits required")
    species=set(r)|set(c)
    missing=sorted(species-set(external_reference))
    if missing:
        raise ValueError("Cannot impute source missing reference trait: "+",".join(missing[:8]))
    ref={}
    for key in species:
        z=external_reference[key]
        if isinstance(z,bool) or not isinstance(z,(int,float)) or not math.isfinite(z):
            raise ValueError("Original reference species trait must be finite")
        ref[key]=float(z)
    base=sum((r.get(k,(0.,0.))[0]-c.get(k,(0.,0.))[0])*ref[k] for k in species)
    within_r=sum(w*(z-ref[k]) for k,(w,z) in r.items())
    within_c=sum(w*(z-ref[k]) for k,(w,z) in c.items())
    rail_cwm=sum(w*z for w,z in r.values())
    control_cwm=sum(w*z for w,z in c.values())
    identity=base+within_r-within_c
    delta=rail_cwm-control_cwm
    if not math.isclose(delta,identity,rel_tol=1.e-10,abs_tol=1.e-10):
        raise ValueError("Exact community trait change accounting identity violated")
    jointly_observed=sorted(set(r)&set(c))
    shifts={k:{
        "source_rail_within_species_trait_mean":r[k][1],
        "source_nonrail_within_species_trait_mean":c[k][1],
        "source_rail_minus_nonrail_within_species_trait_shift":r[k][1]-c[k][1],
        "external_independent_species_reference_trait":ref[k],
    } for k in jointly_observed}
    return {
       "rail_CWM":rail_cwm,"control_CWM":control_cwm,"observed_rail_minus_control_CWM":delta,
       "species_turnover_component_relative_to_reference":base,
       "rail_site_within_species_reference_deviation":within_r,
       "control_site_within_species_reference_deviation":within_c,
       "within_species_component_rail_minus_control":within_r-within_c,
       "identity_reconstruction":identity,
       "n_species_observed_in_rail":len(r),
       "n_species_observed_in_control":len(c),
       "n_jointly_observed_with_real_within_species_comparison":len(jointly_observed),
       "source_jointly_observed_species_trait_shifts":shifts,
       "reference_trait_basis_is_external_and_must_be_provenance_verified":True,
       "weights_original_observational_measure_and_matched_effort_not_checked_by_math":True,
       "plasticity_vs_genetic_selection_NOT_identified_from_CWM":True,
       "no_rail_transport_or_photophenology_effect_inferred":True,
    }
