"""Sealed end-to-end bridge from replayed WLS centres and R2 disks to T1.

The individual T1, R2-component, and margin-transfer certificates deliberately
have narrow scopes.  This module supplies the missing composition contract.  A
bundle is issued only when all of the following replay successfully:

* the complete frozen R2 protocol is a deterministic PASS;
* one replayed digital payload and WLS execution record fixes the centre of
  every canonical rate--harmonic disk of the old T1 problem;
* one certified R2 component-bound request fixes the non-numerical radius,
  while the target-specific WLS solve enclosure is added to its sampling floor;
* units, split identifiers, estimator geometry, Gram evidence, and component
  radii agree across those artifacts;
* the R2 result rows construct a complete new response-disk family; and
* the old robust T1 witness survives that family under the one-way margin
  transfer theorem.

The WLS execution v3 enclosure alone is deliberately narrower than an
exact-arithmetic sample-to-coefficient proof.  Bundle v3 therefore also
requires a replayed exact integer-payload interval certificate covering phase,
``G``/``b`` accumulation, the joint solve, and front-end division.  The full
published complex-disk radius replaces the execution record's solve-only floor;
the two numerical floors are never added twice.

Every failure is ``NOT_CERTIFIABLE``.  This module has no outer or rejection
relation, so a broken evidence chain cannot be misread as model falsification.
"""

from __future__ import annotations

import hashlib
import json
import math
import struct
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, replace
from enum import Enum
from fractions import Fraction
from typing import Any

import numpy as np

from .certified_set import (
    CertifiedHarmonicData,
    CertifiedJointProblem,
    CertifiedParameterBox,
)
from .margin_transfer import (
    DiskMarginPerturbation,
    MarginTransferRelation,
    certify_margin_transfer,
    replay_margin_transfer_certificate,
)
from .r2_bound_certificate import (
    R2BoundCertificate,
    R2BoundRelation,
    R2BoundRequest,
    replay_r2_bound_certificate,
)
from .r2_bound_certificate import (
    _safe_manifest as _r2_request_manifest,
)
from .r2_protocol import (
    R2ProtocolStatus,
    canonical_protocol_sha256,
    canonical_window_result_sha256,
    validate_r2_protocol,
)
from .tube_certificate import TubeRelation, _problem_manifest, replay_tube_certificate
from .wls_execution import (
    SCHEMA_VERSION as WLS_EXECUTION_SCHEMA,
)
from .wls_execution import (
    bounded_sample_indices,
    replay_wls_execution_record,
)
from .wls_interval_certificate import (
    SCHEMA_VERSION as WLS_INTERVAL_SCHEMA,
)
from .wls_interval_certificate import (
    replay_wls_interval_certificate,
)

SCHEMA = "fo-ekf.r2-t1-evidence-bundle.v3"
ALGORITHM = "exact-sample-wls-centre-and-r2-disk-to-frozen-t1-selector-v3"
MAX_CERTIFICATE_BYTES = 32_000_000
MAX_CHILD_CERTIFICATE_BYTES = 8_000_000
MAX_COMPONENT_LINKS = 4096
MAX_PROTOCOL_NODES = 200_000
MAX_PROTOCOL_DEPTH = 64
MAX_PROTOCOL_TEXT_BYTES = 16_000_000
MAX_WLS_SAMPLES = 4096


class R2T1BundleRelation(str, Enum):
    """The only conclusions expressible by the evidence bundle."""

    PRESERVED_ROBUST = "PRESERVED_ROBUST"
    PRESERVED_CLOSED = "PRESERVED_CLOSED"
    NOT_CERTIFIABLE = "NOT_CERTIFIABLE"


@dataclass(frozen=True)
class R2DiskCertificateLink:
    """One replayed WLS centre and component certificate for a T1 disk.

    ``None`` defaults let callers load legacy links safely, but v3 issuance
    rejects them explicitly instead of silently accepting an unbound centre.
    """

    harmonic_index: int
    rate_id: str
    request: R2BoundRequest
    certificate_json: str | R2BoundCertificate
    wls_execution_record: dict[str, Any] | None = None
    digital_samples: tuple[int, ...] | None = None
    wls_interval_certificate: dict[str, Any] | None = None

    def json_text(self) -> str:
        value = self.certificate_json
        return value.to_json() if isinstance(value, R2BoundCertificate) else value


@dataclass(frozen=True)
class ResponseDiskUpdate:
    """A deterministic R2 response disk and its old-centre displacement."""

    harmonic_index: int
    rate_id: str
    window_id: str
    center: complex
    radius: float
    center_shift_abs_upper: float


@dataclass(frozen=True)
class R2T1BundleCertificate:
    """A sealed composition certificate or fail-closed evidence record."""

    relation: R2T1BundleRelation
    reason: str
    input_sha256: str
    certificate_sha256: str
    certificate_json: str
    response_disks: tuple[ResponseDiskUpdate, ...] = ()
    updated_problem: CertifiedJointProblem | None = None

    def to_json(self) -> str:
        return self.certificate_json


@dataclass(frozen=True)
class R2T1BundleReplayResult:
    """Replay verdict for a bundle and its reconstructed response family."""

    valid: bool
    relation: R2T1BundleRelation
    reason: str
    response_disks: tuple[ResponseDiskUpdate, ...] = ()
    updated_problem: CertifiedJointProblem | None = None


@dataclass(frozen=True)
class _Evaluation:
    relation: R2T1BundleRelation
    reason: str
    input_manifest: dict[str, Any]
    proof: dict[str, Any] | None
    response_disks: tuple[ResponseDiskUpdate, ...] = ()
    updated_problem: CertifiedJointProblem | None = None


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    )


def _sha256_json(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _sha256_text(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def _plain_tree_resource_reason(value: Any) -> str | None:
    """Bound untrusted protocol trees without recursive traversal or repr."""

    stack: list[tuple[Any, int]] = [(value, 0)]
    nodes = 0
    text_bytes = 0
    while stack:
        item, depth = stack.pop()
        nodes += 1
        if nodes > MAX_PROTOCOL_NODES:
            return "protocol_resource_nodes"
        if depth > MAX_PROTOCOL_DEPTH:
            return "protocol_resource_depth"
        if type(item) is dict:
            for key, child in item.items():
                if type(key) is not str:
                    return "protocol_key_type_invalid"
                text_bytes += len(key.encode("utf-8"))
                stack.append((child, depth + 1))
        elif type(item) is list:
            stack.extend((child, depth + 1) for child in item)
        elif type(item) is str:
            text_bytes += len(item.encode("utf-8"))
        elif item is None or type(item) in {bool, int, float}:
            pass
        else:
            return "protocol_value_type_invalid"
        if text_bytes > MAX_PROTOCOL_TEXT_BYTES:
            return "protocol_resource_text_bytes"
    return None


def _safe_type_label(value: Any) -> str:
    value_type = type(value)
    return f"{value_type.__module__}.{value_type.__qualname__}"


def _safe_problem_digest(problem: Any, box: Any) -> str:
    if type(problem) is not CertifiedJointProblem or type(box) is not CertifiedParameterBox:
        return _sha256_json(
            {
                "problem_type": _safe_type_label(problem),
                "box_type": _safe_type_label(box),
            }
        )
    try:
        if (
            type(problem.harmonics) is not tuple
            or len(problem.harmonics) > MAX_COMPONENT_LINKS
            or type(box.damping_offsets) is not tuple
            or len(box.damping_offsets) > MAX_COMPONENT_LINKS
        ):
            return "resource_limit"
        return _sha256_json(_problem_manifest(problem, box))
    except (AttributeError, OverflowError, TypeError, ValueError, RecursionError):
        return "invalid"


def _safe_request_digest(request: Any) -> str:
    if type(request) is not R2BoundRequest:
        return _sha256_json({"request_type": _safe_type_label(request)})
    try:
        if (
            type(request.sample_times_s) is not tuple
            or len(request.sample_times_s) > 4096
            or (
                request.sample_indices is not None
                and (
                    type(request.sample_indices) is not tuple or len(request.sample_indices) > 4096
                )
            )
            or type(request.retained_harmonics) is not tuple
            or len(request.retained_harmonics) > 64
        ):
            return "resource_limit"
        return _sha256_json(_r2_request_manifest(request))
    except (AttributeError, OverflowError, TypeError, ValueError, RecursionError):
        return "invalid"


def _rat(value: Fraction) -> str:
    value = Fraction(value)
    return (
        str(value.numerator) if value.denominator == 1 else f"{value.numerator}/{value.denominator}"
    )


def _number_fraction(value: Any) -> Fraction:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError("finite numeric value required")
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("finite numeric value required")
        return Fraction.from_float(value)
    return Fraction(value)


def _same_number(left: Any, right: Any) -> bool:
    try:
        return _number_fraction(left) == _number_fraction(right)
    except ValueError:
        return False


def _parse_proof_fraction(value: Any) -> Fraction:
    if not isinstance(value, str) or not value or len(value) > 20_000:
        raise ValueError("invalid proof rational")
    parsed = Fraction(value)
    if _rat(parsed) != value:
        raise ValueError("noncanonical proof rational")
    return parsed


def _load_json_object(value: str, *, byte_limit: int) -> dict[str, Any]:
    if not isinstance(value, str) or len(value.encode("utf-8")) > byte_limit:
        raise ValueError("JSON evidence exceeds resource limit")

    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, item in pairs:
            if key in result:
                raise ValueError("duplicate JSON key")
            result[key] = item
        return result

    try:
        result = json.loads(value, object_pairs_hook=reject_duplicates)
    except RecursionError as exc:
        raise ValueError("JSON evidence exceeds nesting limit") from exc
    if not isinstance(result, dict):
        raise ValueError("JSON evidence must be an object")
    return result


def _seal(payload: dict[str, Any]) -> tuple[str, str]:
    digest = _sha256_json(payload)
    document = dict(payload)
    document["certificate_sha256"] = digest
    return _canonical_json(document), digest


def _certificate_text(link: Any) -> str:
    if type(link) is not R2DiskCertificateLink:
        raise ValueError("component link type invalid")
    value = link.certificate_json
    if type(value) is R2BoundCertificate:
        value = value.to_json()
    if type(value) is not str:
        raise ValueError("component certificate must be JSON text")
    return value


def _digital_payload_sha256(samples: Any) -> str:
    """Hash an exact little-endian int64 payload without unsafe coercion."""

    if type(samples) is not tuple or not 1 <= len(samples) <= MAX_WLS_SAMPLES:
        raise ValueError("digital payload shape invalid")
    digest = hashlib.sha256()
    for value in samples:
        if type(value) is not int or not -(2**63) <= value < 2**63:
            raise ValueError("digital payload value invalid")
        digest.update(struct.pack("<q", value))
    return digest.hexdigest()


def _safe_wls_record_digest(record: Any) -> tuple[str, str | None]:
    if type(record) is not dict:
        return _sha256_json({"record_type": _safe_type_label(record)}), "wls_record_type_invalid"
    resource_reason = _plain_tree_resource_reason(record)
    if resource_reason is not None:
        wls_reason = f"wls_execution_{resource_reason.removeprefix('protocol_')}"
        return f"invalid:{wls_reason}", wls_reason
    try:
        encoded = _canonical_json(record).encode("utf-8")
    except (OverflowError, TypeError, ValueError, RecursionError):
        return "invalid", "wls_record_json_invalid"
    if len(encoded) > MAX_CHILD_CERTIFICATE_BYTES:
        return "invalid:wls_record_resource_bytes", "wls_record_resource_bytes"
    return hashlib.sha256(encoded).hexdigest(), None


def _safe_interval_certificate_digest(certificate: Any) -> tuple[str, str | None]:
    if type(certificate) is not dict:
        return (
            _sha256_json({"certificate_type": _safe_type_label(certificate)}),
            "wls_interval_certificate_type_invalid",
        )
    resource_reason = _plain_tree_resource_reason(certificate)
    if resource_reason is not None:
        reason = f"wls_interval_{resource_reason.removeprefix('protocol_')}"
        return f"invalid:{reason}", reason
    try:
        encoded = _canonical_json(certificate).encode("utf-8")
    except (OverflowError, TypeError, ValueError, RecursionError):
        return "invalid", "wls_interval_certificate_json_invalid"
    if len(encoded) > MAX_CHILD_CERTIFICATE_BYTES:
        return "invalid:wls_interval_certificate_resource_bytes", (
            "wls_interval_certificate_resource_bytes"
        )
    return hashlib.sha256(encoded).hexdigest(), None


def _safe_manifest(
    problem: Any,
    box: Any,
    protocol_config: Any,
    baseline_certificate_json: Any,
    links: Any,
) -> dict[str, Any]:
    """Create an inert input binding even when issuance inputs are malformed."""

    protocol_resource_reason = (
        _plain_tree_resource_reason(protocol_config)
        if type(protocol_config) is dict
        else "protocol_config_type_invalid"
    )
    try:
        protocol_hash = (
            canonical_protocol_sha256(protocol_config)
            if protocol_resource_reason is None
            else f"invalid:{protocol_resource_reason}"
        )
    except (OverflowError, TypeError, ValueError, RecursionError):
        protocol_hash = "invalid"
    baseline_hash = (
        _sha256_text(baseline_certificate_json)
        if type(baseline_certificate_json) is str
        else _sha256_json({"baseline_type": _safe_type_label(baseline_certificate_json)})
    )
    rows: list[dict[str, Any]] = []
    if type(links) in {list, tuple}:
        for link in links[:MAX_COMPONENT_LINKS]:
            if type(link) is R2DiskCertificateLink:
                try:
                    text = _certificate_text(link)
                    certificate_hash = _sha256_text(text)
                except (TypeError, ValueError):
                    certificate_hash = "invalid"
                wls_record_hash, wls_record_reason = _safe_wls_record_digest(
                    link.wls_execution_record
                )
                try:
                    digital_payload_hash = _digital_payload_sha256(link.digital_samples)
                except (OverflowError, TypeError, ValueError):
                    digital_payload_hash = "invalid"
                interval_hash, interval_reason = _safe_interval_certificate_digest(
                    link.wls_interval_certificate
                )
                declared_exact_coverage = (
                    type(link.wls_interval_certificate) is dict
                    and link.wls_interval_certificate.get("schema_version") == WLS_INTERVAL_SCHEMA
                    and link.wls_interval_certificate.get("relation") == "CERTIFIED_INTERVAL"
                )
                rows.append(
                    {
                        "harmonic_index": link.harmonic_index,
                        "rate_id": link.rate_id,
                        "request_manifest_sha256": _safe_request_digest(link.request),
                        "certificate_json_sha256": certificate_hash,
                        "wls_execution_record_json_sha256": wls_record_hash,
                        "wls_execution_record_resource_reason": wls_record_reason,
                        "external_digital_payload_sha256": digital_payload_hash,
                        "wls_interval_certificate_json_sha256": interval_hash,
                        "wls_interval_certificate_resource_reason": interval_reason,
                        "declared_exact_sample_wls_coverage": declared_exact_coverage,
                    }
                )
            else:
                rows.append({"malformed_link_type": _safe_type_label(link)})
        if len(links) > MAX_COMPONENT_LINKS:
            rows.append({"component_links_truncated_at": MAX_COMPONENT_LINKS})
    else:
        rows.append({"malformed_links_type": _safe_type_label(links)})
    return {
        "semantics": (
            "conditional-declared-r2-disks-with-replayed-binary64-wls-centres;"
            "high-level-preservation-requires-exact-sample-wls-enclosure"
        ),
        "mathematical_exact_sample_wls_coverage": "requires_replay",
        "uncovered_numeric_term": "requires_replay",
        "problem_and_parameter_box_sha256": _safe_problem_digest(problem, box),
        "protocol_sha256": protocol_hash,
        "protocol_resource_reason": protocol_resource_reason,
        "baseline_t1_certificate_json_sha256": baseline_hash,
        "component_links": rows,
    }


def _fail(manifest: dict[str, Any], reason: str) -> _Evaluation:
    return _Evaluation(R2T1BundleRelation.NOT_CERTIFIABLE, reason, manifest, None)


def _evaluation_payload(evaluation: _Evaluation) -> dict[str, Any]:
    payload: dict[str, Any] = {
        "schema": SCHEMA,
        "algorithm": ALGORITHM,
        "input_manifest": evaluation.input_manifest,
        "input_sha256": _sha256_json(evaluation.input_manifest),
        "relation": evaluation.relation.value,
        "reason": evaluation.reason,
    }
    if evaluation.proof is not None:
        payload["proof"] = evaluation.proof
    return payload


def _apply_bundle_size_limit(
    evaluation: _Evaluation,
    fallback_manifest: dict[str, Any],
) -> _Evaluation:
    candidate_json, _ = _seal(_evaluation_payload(evaluation))
    if len(candidate_json.encode("utf-8")) <= MAX_CERTIFICATE_BYTES:
        return evaluation
    return _fail(fallback_manifest, "bundle_certificate_resource_bytes")


def _exact_hypot_upper(real: Fraction, imag: Fraction) -> Fraction:
    """Return an exact binary rational proven above ``sqrt(real^2+imag^2)``."""

    square = real * real + imag * imag
    if square == 0:
        return Fraction(0)
    candidate = math.nextafter(math.sqrt(float(square)), math.inf)
    if not math.isfinite(candidate):
        raise ValueError("centre displacement is nonfinite")
    for _ in range(64):
        bound = Fraction.from_float(candidate)
        if bound >= 0 and bound * bound >= square:
            return bound
        candidate = math.nextafter(candidate, math.inf)
    raise ValueError("centre displacement could not be bounded")


def _expected_layout(problem: CertifiedJointProblem) -> tuple[tuple[int, str], ...]:
    return tuple(
        (harmonic.harmonic_index, rate_id)
        for harmonic in problem.harmonics
        for rate_id in problem.rate_ids
    )


def _protocol_table(config: Mapping[str, Any], name: str) -> Mapping[str, Any]:
    value = config.get(name)
    if not isinstance(value, Mapping):
        raise ValueError(f"protocol table {name} missing")
    return value


def _baseline_model_contract_reason(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    config: Mapping[str, Any],
) -> str | None:
    model = _protocol_table(config, "model")
    estimator = _protocol_table(config, "estimator")
    if model.get("t1_rate_ids") != list(problem.rate_ids):
        return "baseline_t1_rate_layout_mismatch"
    if model.get("damping_anchor_rate_id") != problem.anchor_rate_id:
        return "baseline_t1_damping_anchor_mismatch"
    alpha = model.get("alpha_interval")
    damping = model.get("lambda_interval")
    if not isinstance(alpha, list) or len(alpha) != 2:
        return "protocol_alpha_interval_invalid"
    if not isinstance(damping, list) or len(damping) != 2:
        return "protocol_lambda_interval_invalid"
    try:
        if not (
            _number_fraction(alpha[0]) <= Fraction.from_float(box.order.lower)
            and Fraction.from_float(box.order.upper) <= _number_fraction(alpha[1])
        ):
            return "baseline_t1_alpha_outside_protocol_model"
        if not (
            _number_fraction(damping[0]) <= Fraction.from_float(box.damping.lower)
            and Fraction.from_float(box.damping.upper) <= _number_fraction(damping[1])
        ):
            return "baseline_t1_damping_outside_protocol_model"
    except ValueError:
        return "baseline_t1_model_interval_invalid"

    retained = estimator.get("retained_harmonics")
    q_lower = model.get("q_abs_lower_by_harmonic")
    q_upper = model.get("q_abs_upper_by_harmonic")
    if not all(isinstance(value, list) for value in (retained, q_lower, q_upper)):
        return "protocol_q_prior_layout_invalid"
    if not (len(retained) == len(q_lower) == len(q_upper)):
        return "protocol_q_prior_layout_invalid"
    for harmonic in problem.harmonics:
        try:
            position = retained.index(harmonic.harmonic_index)
        except ValueError:
            return "baseline_t1_harmonic_absent_from_protocol"
        if not _same_number(q_lower[position], harmonic.morphology_bounds[0]):
            return "baseline_t1_q_lower_mismatch"
        if not _same_number(q_upper[position], harmonic.morphology_bounds[1]):
            return "baseline_t1_q_upper_mismatch"
    return None


def _component_contract_reason(
    *,
    request: R2BoundRequest,
    row: Mapping[str, Any],
    config: Mapping[str, Any],
    protocol_hash: str,
    expected_harmonic: int,
    expected_rate: str,
    expected_base_frequency: float,
) -> str | None:
    data = _protocol_table(config, "data_boundary")
    units = _protocol_table(config, "units")
    model = _protocol_table(config, "model")
    estimator = _protocol_table(config, "estimator")
    preprocessing = _protocol_table(config, "preprocessing")
    bounds = _protocol_table(config, "bounds")
    sampling_table = bounds.get("sampling")
    if not isinstance(sampling_table, Mapping):
        return "protocol_sampling_table_missing"

    if request.protocol_sha256 != protocol_hash:
        return "component_protocol_hash_mismatch"
    if request.frozen_before_target_ecg_access is not True:
        return "component_not_frozen_before_target_access"
    if request.signal_unit != units.get("signal_unit"):
        return "component_signal_unit_mismatch"
    if request.calibration_split_id != data.get("calibration_split_id"):
        return "component_calibration_split_mismatch"
    if request.identification_split_id != data.get("identification_split_id"):
        return "component_identification_split_mismatch"
    if request.calibration_split_id == request.identification_split_id:
        return "component_split_overlap"
    if not _same_number(request.tau_star_s, units.get("tau_star_s")):
        return "component_tau_star_mismatch"
    if list(request.alpha_interval) != model.get("alpha_interval"):
        return "component_alpha_interval_mismatch"
    if list(request.lambda_interval) != model.get("lambda_interval"):
        return "component_lambda_interval_mismatch"
    if list(request.retained_harmonics) != estimator.get("retained_harmonics"):
        return "component_harmonic_basis_mismatch"
    if not _same_number(
        request.gram_min_eigenvalue_threshold,
        estimator.get("gram_min_eigenvalue_threshold"),
    ):
        return "component_gram_min_threshold_mismatch"
    if not _same_number(
        request.gram_condition_number_max,
        estimator.get("gram_condition_number_max"),
    ):
        return "component_gram_condition_threshold_mismatch"
    expected_front_end = (
        "identity"
        if preprocessing.get("identity_front_end") is True
        else "preconvolved_output_envelopes"
    )
    if request.dynamic_front_end_mode != expected_front_end:
        return "component_front_end_semantics_mismatch"
    if request.sampling is None or request.sampling.implementation != sampling_table.get(
        "implementation"
    ):
        return "component_estimator_implementation_mismatch"
    if request.sampling.quadrature_decomposition_target != sampling_table.get(
        "quadrature_decomposition_target"
    ):
        return "component_quadrature_decomposition_target_mismatch"

    if row.get("record_manifest_sha256") != request.record_manifest_sha256:
        return "component_record_manifest_hash_mismatch"
    if row.get("estimator_manifest_sha256") != request.estimator_manifest_sha256:
        return "component_estimator_manifest_hash_mismatch"
    evidence_groups = {
        "history": request.history,
        "dwell": request.dwell,
        "sampling": request.sampling,
        "delay": request.delay,
    }
    for name, component_group in evidence_groups.items():
        protocol_group = bounds.get(name)
        if component_group is None or not isinstance(protocol_group, Mapping):
            return f"component_{name}_evidence_missing"
        evidence = component_group.evidence
        for protocol_key, component_value in (
            ("source_kind", evidence.source_kind),
            ("source_reference", evidence.source_reference),
            ("calibration_split_id", evidence.calibration_split_id),
            ("independent_of_identification_fit", evidence.independent_of_identification_fit),
            ("uses_target_fit_residuals", evidence.uses_target_fit_residuals),
            ("evidence_sha256", evidence.evidence_sha256),
        ):
            if protocol_group.get(protocol_key) != component_value:
                return f"component_{name}_{protocol_key}_mismatch"

    if request.target_harmonic != expected_harmonic:
        return "component_target_harmonic_mismatch"
    if request.window_id != row.get("window_id"):
        return "component_window_id_mismatch"
    if row.get("rate_id") != expected_rate or row.get("target_harmonic") != expected_harmonic:
        return "window_disk_identity_mismatch"
    if not _same_number(row.get("nominal_rate_rad_s"), expected_base_frequency):
        return "t1_r2_frequency_mismatch"
    for request_value, row_key in (
        (request.left_r_peak_time_s, "left_r_peak_time_s"),
        (request.right_r_peak_time_s, "right_r_peak_time_s"),
        (request.complete_rr_intervals, "complete_rr_intervals"),
    ):
        if not _same_number(request_value, row.get(row_key)):
            return f"component_{row_key}_mismatch"
    if row.get("sample_count") != len(request.sample_times_s):
        return "component_sample_count_mismatch"
    return None


def _float64_payload_sha256(values: Sequence[Any]) -> str:
    digest = hashlib.sha256()
    for value in values:
        fraction = _number_fraction(value)
        encoded = float(value)
        if not math.isfinite(encoded) or Fraction.from_float(encoded) != fraction:
            raise ValueError("value is not an exact finite binary64 number")
        digest.update(struct.pack("<d", encoded))
    return digest.hexdigest()


def _parse_wls_dyadic_upper(value: Any) -> Fraction:
    if type(value) is not dict:
        raise ValueError("WLS dyadic bound must be a table")
    numerator_text = value.get("dyadic_numerator")
    denominator_text = value.get("dyadic_denominator")
    float_upper = value.get("float_upper")
    if type(numerator_text) is not str or type(denominator_text) is not str:
        raise ValueError("WLS dyadic integers must be strings")
    if not numerator_text or not denominator_text or len(numerator_text) > 20_000:
        raise ValueError("WLS dyadic integer invalid")
    numerator = int(numerator_text)
    denominator = int(denominator_text)
    if str(numerator) != numerator_text or str(denominator) != denominator_text or denominator <= 0:
        raise ValueError("WLS dyadic integer is noncanonical")
    bound = Fraction(numerator, denominator)
    if bound < 0 or _number_fraction(float_upper) != bound:
        raise ValueError("WLS float and dyadic upper bounds disagree")
    return bound


def _parse_interval_dyadic_upper(value: Any) -> Fraction:
    if type(value) is not dict or type(value.get("dyadic")) is not dict:
        raise ValueError("interval upper bound must contain a dyadic table")
    dyadic = value["dyadic"]
    numerator_text = dyadic.get("numerator")
    denominator_text = dyadic.get("denominator")
    if type(numerator_text) is not str or type(denominator_text) is not str:
        raise ValueError("interval dyadic integers must be strings")
    numerator = int(numerator_text)
    denominator = int(denominator_text)
    if str(numerator) != numerator_text or str(denominator) != denominator_text or denominator <= 0:
        raise ValueError("interval dyadic integer is noncanonical")
    bound = Fraction(numerator, denominator)
    if bound < 0 or _number_fraction(value.get("float_upper")) != bound:
        raise ValueError("interval float and dyadic upper bounds disagree")
    return bound


def _wls_execution_contract(
    *,
    link: R2DiskCertificateLink,
    request: R2BoundRequest,
    row: Mapping[str, Any],
    config: Mapping[str, Any],
    protocol_hash: str,
) -> tuple[str | None, Fraction, dict[str, Any] | None]:
    """Replay and cross-bind one external digital payload to one R2 centre."""

    record = link.wls_execution_record
    samples = link.digital_samples
    if record is None:
        return "wls_execution_record_missing", Fraction(0), None
    if samples is None:
        return "wls_digital_payload_missing", Fraction(0), None
    record_digest, record_resource_reason = _safe_wls_record_digest(record)
    if record_resource_reason is not None:
        return record_resource_reason, Fraction(0), None
    try:
        digital_digest = _digital_payload_sha256(samples)
    except (OverflowError, TypeError, ValueError):
        return "wls_digital_payload_invalid", Fraction(0), None
    # A mismatch may be rejected from the inert, resource-bounded record before
    # the comparatively expensive numerical replay.  These checks never grant
    # coverage; the matching path still requires full execution replay below.
    try:
        preflight_window = record["window"]
        if type(preflight_window) is not dict:
            raise ValueError("WLS window table missing")
        preflight_start = preflight_window.get("start_sample_inclusive")
        preflight_end = preflight_window.get("end_sample_exclusive")
        if (
            request.window_start_sample != preflight_start
            or request.window_end_sample_exclusive != preflight_end
            or row.get("wls_window_start_sample") != preflight_start
            or row.get("wls_window_end_sample_exclusive") != preflight_end
        ):
            return "wls_r2_integer_window_geometry_mismatch", Fraction(0), None
        preflight_indices = tuple(
            int(value) for value in bounded_sample_indices(preflight_start, preflight_end)
        )
        if request.sample_indices != preflight_indices:
            return "wls_r2_exact_sample_indices_mismatch", Fraction(0), None
    except (IndexError, KeyError, OverflowError, TypeError, ValueError):
        return "wls_integer_phase_preflight_invalid", Fraction(0), None
    try:
        replay = replay_wls_execution_record(record, samples)
    except (
        ArithmeticError,
        IndexError,
        KeyError,
        OverflowError,
        TypeError,
        ValueError,
        RecursionError,
    ):
        return "wls_execution_replay_failed", Fraction(0), None
    if replay.get("status") != "MATCH":
        return "wls_execution_record_not_replayable", Fraction(0), None

    try:
        if record.get("schema_version") != WLS_EXECUTION_SCHEMA:
            return "wls_execution_schema_unsupported", Fraction(0), None
        binding = record["binding"]
        window = record["window"]
        payload_hashes = record["payload_hashes"]
        weights = record["weights"]
        front_end = record["front_end"]
        result = record["result"]
        if not all(
            type(value) is dict
            for value in (binding, window, payload_hashes, weights, front_end, result)
        ):
            raise ValueError("WLS record table missing")

        if binding.get("protocol_sha256") != protocol_hash:
            return "wls_protocol_hash_mismatch", Fraction(0), None
        if binding.get("record_sha256") != request.record_manifest_sha256:
            return "wls_record_manifest_hash_mismatch", Fraction(0), None
        if binding.get("estimator_sha256") != request.estimator_manifest_sha256:
            return "wls_estimator_manifest_hash_mismatch", Fraction(0), None
        if row.get("record_manifest_sha256") != binding.get("record_sha256"):
            return "window_wls_record_hash_mismatch", Fraction(0), None
        if row.get("estimator_manifest_sha256") != binding.get("estimator_sha256"):
            return "window_wls_estimator_hash_mismatch", Fraction(0), None
        if row.get("wls_window_sha256") != binding.get("window_sha256"):
            return "window_wls_window_hash_mismatch", Fraction(0), None
        if row.get("wls_execution_sha256") != record.get("execution_sha256"):
            return "window_wls_execution_hash_mismatch", Fraction(0), None

        expected_payload_rows = {
            "wls_sample_index_sha256": "sample_index_sha256",
            "wls_digital_sample_sha256": "digital_sample_sha256",
            "wls_time_seconds_sha256": "time_seconds_sha256",
            "wls_normalized_weights_sha256": "normalized_weights_sha256",
        }
        for row_key, record_key in expected_payload_rows.items():
            if row.get(row_key) != payload_hashes.get(record_key):
                return f"window_{row_key}_mismatch", Fraction(0), None
        if payload_hashes.get("digital_sample_sha256") != digital_digest:
            return "wls_external_digital_payload_hash_mismatch", Fraction(0), None

        sample_count = window.get("sample_count")
        if type(sample_count) is not int or sample_count != len(samples):
            return "wls_sample_count_mismatch", Fraction(0), None
        if sample_count != len(request.sample_times_s) or sample_count != row.get("sample_count"):
            return "wls_r2_sample_count_mismatch", Fraction(0), None
        start_sample = window.get("start_sample_inclusive")
        end_sample = window.get("end_sample_exclusive")
        if (
            request.window_start_sample != start_sample
            or request.window_end_sample_exclusive != end_sample
            or row.get("wls_window_start_sample") != start_sample
            or row.get("wls_window_end_sample_exclusive") != end_sample
        ):
            return "wls_r2_integer_window_geometry_mismatch", Fraction(0), None
        expected_indices = tuple(
            int(value) for value in bounded_sample_indices(start_sample, end_sample)
        )
        if request.sample_indices != expected_indices:
            return "wls_r2_exact_sample_indices_mismatch", Fraction(0), None
        if payload_hashes.get("sample_index_sha256") != _digital_payload_sha256(
            request.sample_indices
        ):
            return "wls_r2_sample_index_payload_mismatch", Fraction(0), None
        sampling_frequency = window.get("sampling_frequency_hz")
        sampling_table = _protocol_table(_protocol_table(config, "bounds"), "sampling")
        if not _same_number(sampling_frequency, sampling_table.get("sample_rate_hz")):
            return "wls_sampling_frequency_mismatch", Fraction(0), None
        if (
            not _same_number(
                float(window["start_sample_inclusive"]) / float(sampling_frequency),
                request.left_r_peak_time_s,
            )
            or not _same_number(
                float(window["end_sample_exclusive"]) / float(sampling_frequency),
                request.right_r_peak_time_s,
            )
            or window.get("rr_interval_count") != request.complete_rr_intervals
        ):
            return "wls_window_geometry_mismatch", Fraction(0), None
        if payload_hashes.get("time_seconds_sha256") != _float64_payload_sha256(
            request.sample_times_s
        ):
            return "wls_r2_sample_times_mismatch", Fraction(0), None

        raw_weights = request.raw_weights
        if (
            type(raw_weights) is not tuple
            or len(raw_weights) != sample_count
            or any(not _same_number(value, 1.0) for value in raw_weights)
        ):
            return "wls_r2_raw_weights_not_identity", Fraction(0), None
        if payload_hashes.get("raw_weights_sha256") != _float64_payload_sha256(raw_weights):
            return "wls_r2_raw_weights_mismatch", Fraction(0), None
        normalized = tuple(1.0 / sample_count for _ in range(sample_count))
        if payload_hashes.get("normalized_weights_sha256") != _float64_payload_sha256(normalized):
            return "wls_r2_normalized_weights_mismatch", Fraction(0), None
        if weights.get("count") != sample_count:
            return "wls_weight_count_mismatch", Fraction(0), None

        if window.get("retained_harmonics") != list(request.retained_harmonics):
            return "wls_retained_harmonics_mismatch", Fraction(0), None
        positive_harmonics = window.get("positive_harmonics")
        expected_positive = [value for value in request.retained_harmonics if value > 0]
        if (
            positive_harmonics != expected_positive
            or request.target_harmonic not in expected_positive
        ):
            return "wls_positive_harmonics_mismatch", Fraction(0), None
        target_position = expected_positive.index(request.target_harmonic)

        h_rows = front_end.get("H_hat_by_harmonic")
        z_rows = result.get("z_tilde")
        if type(h_rows) is not list or type(z_rows) is not list:
            raise ValueError("WLS target rows missing")
        if [value.get("harmonic") for value in h_rows if type(value) is dict] != expected_positive:
            return "wls_front_end_harmonic_layout_mismatch", Fraction(0), None
        if [value.get("harmonic") for value in z_rows if type(value) is dict] != expected_positive:
            return "wls_target_harmonic_layout_mismatch", Fraction(0), None
        h_target = h_rows[target_position]
        z_target = z_rows[target_position]
        if not _same_number(
            row.get("nominal_transfer_real"), h_target.get("real")
        ) or not _same_number(row.get("nominal_transfer_imag"), h_target.get("imag")):
            return "wls_front_end_target_mismatch", Fraction(0), None
        expected_mode = (
            "identity" if request.dynamic_front_end_mode == "identity" else "caller_supplied"
        )
        if front_end.get("mode") != expected_mode:
            return "wls_front_end_mode_mismatch", Fraction(0), None
        if not _same_number(row.get("z_tilde_real"), z_target.get("real")) or not _same_number(
            row.get("z_tilde_imag"), z_target.get("imag")
        ):
            return "wls_target_center_mismatch", Fraction(0), None

    except (IndexError, KeyError, OverflowError, TypeError, ValueError, ZeroDivisionError):
        return "wls_execution_contract_invalid", Fraction(0), None

    interval_certificate = link.wls_interval_certificate
    if interval_certificate is None:
        return "wls_interval_certificate_missing", Fraction(0), None
    interval_digest, interval_resource_reason = _safe_interval_certificate_digest(
        interval_certificate
    )
    if interval_resource_reason is not None:
        return interval_resource_reason, Fraction(0), None
    try:
        digital_array = np.asarray(samples, dtype=np.int64)
        interval_replay = replay_wls_interval_certificate(
            interval_certificate,
            record,
            digital_array,
        )
    except (
        ArithmeticError,
        IndexError,
        KeyError,
        OverflowError,
        TypeError,
        ValueError,
        RecursionError,
    ):
        return "wls_interval_certificate_replay_failed", Fraction(0), None
    if interval_replay.get("status") != "MATCH":
        return "wls_interval_certificate_not_replayable", Fraction(0), None
    if interval_replay.get("relation") != "CERTIFIED_INTERVAL":
        return "wls_interval_certificate_not_certified", Fraction(0), None

    try:
        if interval_certificate.get("schema_version") != WLS_INTERVAL_SCHEMA:
            return "wls_interval_certificate_schema_unsupported", Fraction(0), None
        interval_binding = interval_certificate["binding"]
        mathematical_contract = interval_certificate["mathematical_contract"]
        interval_disks = interval_certificate["response_disks"]
        if (
            type(interval_binding) is not dict
            or type(mathematical_contract) is not dict
            or type(interval_disks) is not list
        ):
            raise ValueError("interval certificate table missing")
        expected_interval_binding = {
            "execution_sha256": record["execution_sha256"],
            "protocol_sha256": protocol_hash,
            "record_sha256": request.record_manifest_sha256,
            "estimator_sha256": request.estimator_manifest_sha256,
            "window_sha256": binding["window_sha256"],
            "digital_sample_sha256": digital_digest,
        }
        if interval_binding != expected_interval_binding:
            return "wls_interval_binding_mismatch", Fraction(0), None
        if row.get("wls_interval_certificate_sha256") != interval_certificate.get(
            "certificate_sha256"
        ):
            return "window_wls_interval_certificate_hash_mismatch", Fraction(0), None
        if row.get("wls_exact_sample_functional_coverage") is not True:
            return "window_wls_exact_coverage_not_declared", Fraction(0), None
        if mathematical_contract.get("positive_harmonics") != expected_positive:
            return "wls_interval_positive_harmonics_mismatch", Fraction(0), None
        if mathematical_contract.get("retained_harmonics") != list(request.retained_harmonics):
            return "wls_interval_retained_harmonics_mismatch", Fraction(0), None
        for field, expected in (
            ("sample_count", sample_count),
            ("start_sample_inclusive", window["start_sample_inclusive"]),
            ("end_sample_exclusive", window["end_sample_exclusive"]),
            ("rr_interval_count", request.complete_rr_intervals),
        ):
            if mathematical_contract.get(field) != expected:
                return f"wls_interval_{field}_mismatch", Fraction(0), None
        if [value.get("harmonic") for value in interval_disks if type(value) is dict] != (
            expected_positive
        ):
            return "wls_interval_response_layout_mismatch", Fraction(0), None
        interval_target = interval_disks[target_position]
        interval_center = interval_target.get("center_binary64")
        if type(interval_center) is not dict:
            raise ValueError("interval target center missing")
        if not _same_number(interval_center.get("real"), z_target.get("real")) or not _same_number(
            interval_center.get("imag"), z_target.get("imag")
        ):
            return "wls_interval_target_center_mismatch", Fraction(0), None
        exact_functional_floor = _parse_interval_dyadic_upper(
            interval_target.get("absolute_error_upper")
        )
        if not _same_number(
            row.get("wls_exact_functional_radius_upper"),
            interval_target["absolute_error_upper"]["float_upper"],
        ):
            return "window_wls_exact_functional_radius_mismatch", Fraction(0), None
        required_coverage = {
            "digital_to_physical_binary64_conversion_rounding",
            "phase_and_complex_exponential_rounding",
            "G_and_b_direct_accumulation_rounding",
            "verified_joint_linear_solve_rounding",
            "front_end_binary64_division_rounding",
        }
        if not required_coverage.issubset(set(interval_certificate.get("coverage", []))):
            return "wls_interval_numeric_coverage_incomplete", Fraction(0), None
        if "G_and_b_accumulation_rounding" in interval_certificate.get("exclusions", []):
            return "wls_interval_accumulation_still_excluded", Fraction(0), None
        if interval_certificate.get("claim_scope") != "frozen_mathematical_sample_functional_only":
            return "wls_interval_claim_scope_mismatch", Fraction(0), None
    except (IndexError, KeyError, OverflowError, TypeError, ValueError, ZeroDivisionError):
        return "wls_interval_certificate_contract_invalid", Fraction(0), None

    evidence = {
        "execution_sha256": record["execution_sha256"],
        "execution_record_json_sha256": record_digest,
        "window_sha256": binding["window_sha256"],
        "external_digital_payload_sha256": digital_digest,
        "sample_index_sha256": payload_hashes["sample_index_sha256"],
        "time_seconds_sha256": payload_hashes["time_seconds_sha256"],
        "normalized_weights_sha256": payload_hashes["normalized_weights_sha256"],
        "target_harmonic": request.target_harmonic,
        "wls_interval_certificate_sha256": interval_certificate["certificate_sha256"],
        "wls_interval_certificate_json_sha256": interval_digest,
        "target_exact_functional_radius_floor": _rat(exact_functional_floor),
        "exact_functional_scope": interval_certificate["claim_scope"],
        "mathematical_exact_sample_wls_coverage": True,
        "uncovered_numeric_term": None,
    }
    return None, exact_functional_floor, evidence


def _component_proof_reason(
    row: Mapping[str, Any],
    component_document: Mapping[str, Any],
    exact_functional_floor: Fraction,
) -> str | None:
    try:
        proof = component_document["proof"]
        gram = proof["gram"]
        dynamic = proof["dynamic"]
        sampling = proof["sampling"]
        delay = proof["delay"]
        digest = component_document["certificate_sha256"]
        if row.get("gram_certificate_sha256") != digest:
            return "window_gram_certificate_hash_mismatch"
        if _number_fraction(row.get("gram_min_eigenvalue")) > _parse_proof_fraction(
            gram["lambda_min_lower"]
        ):
            return "window_gram_min_overclaims_component_certificate"
        if _number_fraction(row.get("gram_condition_number")) < _parse_proof_fraction(
            gram["condition_number_upper"]
        ):
            return "window_gram_condition_understates_component_certificate"
        if _number_fraction(row.get("wls_row_norm_upper_bound")) < _parse_proof_fraction(
            gram["wls_row_norm_upper"]
        ):
            return "window_wls_norm_understates_component_certificate"

        component_floors = {
            "radius_history": _parse_proof_fraction(dynamic["history"]["coefficient_radius_upper"]),
            "radius_dwell": _parse_proof_fraction(dynamic["dwell"]["coefficient_radius_upper"]),
            "radius_sampling": _parse_proof_fraction(sampling["coefficient_radius_upper"]),
            "radius_quadrature": _parse_proof_fraction(
                sampling["quadrature_coefficient_radius_upper"]
            ),
            "radius_delay": _parse_proof_fraction(delay["coefficient_radius_upper"]),
        }
        row_floors = dict(component_floors)
        row_floors["radius_sampling"] += exact_functional_floor
        for field, floor in row_floors.items():
            if _number_fraction(row.get(field)) < floor:
                return f"window_{field}_below_component_certificate"
        if _parse_proof_fraction(proof["partial_radius_upper"]) != sum(
            component_floors.values(), Fraction(0)
        ):
            return "component_partial_radius_sum_mismatch"
    except (KeyError, TypeError, ValueError, ZeroDivisionError):
        return "component_proof_material_invalid"
    return None


def _updated_problem(
    problem: CertifiedJointProblem,
    updates: tuple[ResponseDiskUpdate, ...],
) -> CertifiedJointProblem:
    by_id = {(value.harmonic_index, value.rate_id): value for value in updates}
    harmonics: list[CertifiedHarmonicData] = []
    for harmonic in problem.harmonics:
        replacements = tuple(
            by_id[(harmonic.harmonic_index, rate_id)] for rate_id in problem.rate_ids
        )
        harmonics.append(
            replace(
                harmonic,
                measured_responses=tuple(value.center for value in replacements),
                response_radii=tuple(value.radius for value in replacements),
            )
        )
    return CertifiedJointProblem(tuple(harmonics), problem.anchor_rate_id)


def _update_payload(
    update: ResponseDiskUpdate,
    old_center: complex,
    old_radius: float,
    morphology_radius: float,
) -> dict[str, Any]:
    return {
        "harmonic_index": update.harmonic_index,
        "rate_id": update.rate_id,
        "window_id": update.window_id,
        "old_center": [
            _rat(Fraction.from_float(old_center.real)),
            _rat(Fraction.from_float(old_center.imag)),
        ],
        "new_center": [
            _rat(Fraction.from_float(update.center.real)),
            _rat(Fraction.from_float(update.center.imag)),
        ],
        "old_response_radius": _rat(Fraction.from_float(old_radius)),
        "new_response_radius": _rat(Fraction.from_float(update.radius)),
        "center_shift_abs_upper": _rat(Fraction.from_float(update.center_shift_abs_upper)),
        "old_morphology_radius": _rat(Fraction.from_float(morphology_radius)),
        "new_morphology_radius": _rat(Fraction.from_float(morphology_radius)),
    }


def _evaluate(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    protocol_config: Mapping[str, Any],
    baseline_certificate_json: str,
    links: Sequence[R2DiskCertificateLink],
) -> _Evaluation:
    safe_manifest = _safe_manifest(
        problem,
        box,
        protocol_config,
        baseline_certificate_json,
        links,
    )
    if type(problem) is not CertifiedJointProblem or type(box) is not CertifiedParameterBox:
        return _fail(safe_manifest, "problem_or_parameter_box_type_invalid")
    if type(protocol_config) is not dict:
        return _fail(safe_manifest, "protocol_config_type_invalid")
    protocol_resource_reason = _plain_tree_resource_reason(protocol_config)
    if protocol_resource_reason is not None:
        return _fail(safe_manifest, protocol_resource_reason)
    if type(baseline_certificate_json) is not str:
        return _fail(safe_manifest, "baseline_certificate_type_invalid")
    if type(links) not in {list, tuple}:
        return _fail(safe_manifest, "component_links_type_invalid")
    if len(links) > MAX_COMPONENT_LINKS:
        return _fail(safe_manifest, "component_links_resource_limit")

    try:
        protocol_validation = validate_r2_protocol(protocol_config)
    except (OverflowError, TypeError, ValueError, RecursionError):
        return _fail(safe_manifest, "protocol_validation_failed")
    if protocol_validation.status is not R2ProtocolStatus.PASS_DETERMINISTIC:
        return _fail(safe_manifest, "protocol_not_deterministic_pass")
    protocol_hash = canonical_protocol_sha256(protocol_config)
    if protocol_validation.canonical_sha256 != protocol_hash:
        return _fail(safe_manifest, "protocol_hash_internal_mismatch")

    baseline_replay = replay_tube_certificate(problem, box, baseline_certificate_json)
    if not baseline_replay.valid:
        return _fail(safe_manifest, "baseline_t1_certificate_invalid")
    if baseline_replay.relation is not TubeRelation.ROBUST_INNER:
        return _fail(safe_manifest, "baseline_t1_not_robust")
    try:
        baseline_document = _load_json_object(
            baseline_certificate_json,
            byte_limit=MAX_CHILD_CERTIFICATE_BYTES,
        )
        baseline_digest = baseline_document["certificate_sha256"]
        baseline_input_hash = baseline_document["input_sha256"]
        protocol_table = _protocol_table(protocol_config, "protocol")
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return _fail(safe_manifest, "baseline_t1_envelope_invalid")
    if protocol_table.get("baseline_t1_certificate_sha256") != baseline_digest:
        return _fail(safe_manifest, "baseline_t1_not_bound_by_frozen_protocol")
    if protocol_table.get("baseline_t1_input_sha256") != baseline_input_hash:
        return _fail(safe_manifest, "baseline_t1_input_not_bound_by_frozen_protocol")
    try:
        model_reason = _baseline_model_contract_reason(problem, box, protocol_config)
    except (KeyError, TypeError, ValueError):
        return _fail(safe_manifest, "baseline_t1_model_contract_invalid")
    if model_reason is not None:
        return _fail(safe_manifest, model_reason)

    expected_layout = _expected_layout(problem)
    frozen_links = tuple(links)
    actual_layout: list[tuple[int, str]] = []
    for link in frozen_links:
        if type(link) is not R2DiskCertificateLink:
            return _fail(safe_manifest, "component_link_type_invalid")
        if (
            isinstance(link.harmonic_index, bool)
            or not isinstance(link.harmonic_index, int)
            or not isinstance(link.rate_id, str)
        ):
            return _fail(safe_manifest, "component_link_identity_invalid")
        actual_layout.append((link.harmonic_index, link.rate_id))
    if tuple(actual_layout) != expected_layout:
        return _fail(safe_manifest, "component_disk_coverage_mismatch")

    result_rows = protocol_config.get("window_result")
    if not isinstance(result_rows, list):
        return _fail(safe_manifest, "protocol_window_results_missing")
    rows_by_key: dict[tuple[str, int], Mapping[str, Any]] = {}
    for row in result_rows:
        if (
            not isinstance(row, Mapping)
            or not isinstance(row.get("window_id"), str)
            or not isinstance(row.get("target_harmonic"), int)
            or isinstance(row.get("target_harmonic"), bool)
        ):
            return _fail(safe_manifest, "protocol_window_result_invalid")
        result_key = (row["window_id"], row["target_harmonic"])
        if result_key in rows_by_key:
            return _fail(safe_manifest, "protocol_window_harmonic_not_unique")
        rows_by_key[result_key] = row

    selected_manifest: list[dict[str, Any]] = []
    embedded_components: list[dict[str, Any]] = []
    replayed_wls_evidence: list[dict[str, Any]] = []
    updates: list[ResponseDiskUpdate] = []
    perturbations: list[DiskMarginPerturbation] = []
    update_proof: list[dict[str, Any]] = []
    position = 0
    for harmonic in problem.harmonics:
        for rate_index, rate_id in enumerate(problem.rate_ids):
            link = frozen_links[position]
            position += 1
            if type(link.request) is not R2BoundRequest:
                return _fail(safe_manifest, "component_request_type_invalid")
            row = rows_by_key.get((link.request.window_id, harmonic.harmonic_index))
            if row is None:
                return _fail(safe_manifest, "component_window_result_missing")
            contract_reason = _component_contract_reason(
                request=link.request,
                row=row,
                config=protocol_config,
                protocol_hash=protocol_hash,
                expected_harmonic=harmonic.harmonic_index,
                expected_rate=rate_id,
                expected_base_frequency=harmonic.base_frequencies[rate_index],
            )
            if contract_reason is not None:
                return _fail(safe_manifest, contract_reason)

            try:
                component_json = _certificate_text(link)
            except (TypeError, ValueError):
                return _fail(safe_manifest, "component_certificate_type_invalid")
            component_replay = replay_r2_bound_certificate(link.request, component_json)
            if not component_replay.valid:
                return _fail(safe_manifest, "component_certificate_invalid")
            if component_replay.relation is not R2BoundRelation.CERTIFIED_BOUND:
                return _fail(safe_manifest, "component_certificate_unknown")
            try:
                component_document = _load_json_object(
                    component_json,
                    byte_limit=MAX_CHILD_CERTIFICATE_BYTES,
                )
            except (TypeError, ValueError, json.JSONDecodeError):
                return _fail(safe_manifest, "component_certificate_envelope_invalid")
            wls_reason, exact_functional_floor, wls_evidence = _wls_execution_contract(
                link=link,
                request=link.request,
                row=row,
                config=protocol_config,
                protocol_hash=protocol_hash,
            )
            if wls_reason is not None or wls_evidence is None:
                return _fail(safe_manifest, wls_reason or "wls_execution_contract_invalid")
            proof_reason = _component_proof_reason(
                row,
                component_document,
                exact_functional_floor,
            )
            if proof_reason is not None:
                return _fail(safe_manifest, proof_reason)

            try:
                center_real = float(row["z_tilde_real"])
                center_imag = float(row["z_tilde_imag"])
                radius = float(row["radius_total_deterministic"])
                if not all(math.isfinite(value) for value in (center_real, center_imag, radius)):
                    raise ValueError("nonfinite response disk")
                if radius < 0:
                    raise ValueError("negative response radius")
                old_center = harmonic.measured_responses[rate_index]
                shift = _exact_hypot_upper(
                    Fraction.from_float(center_real) - Fraction.from_float(old_center.real),
                    Fraction.from_float(center_imag) - Fraction.from_float(old_center.imag),
                )
                shift_float = float(shift)
                if Fraction.from_float(shift_float) != shift:
                    raise ValueError("displacement bound lost exactness")
            except (KeyError, OverflowError, TypeError, ValueError):
                return _fail(safe_manifest, "response_disk_construction_failed")

            update = ResponseDiskUpdate(
                harmonic.harmonic_index,
                rate_id,
                link.request.window_id,
                complex(center_real, center_imag),
                radius,
                shift_float,
            )
            morphology_radius = harmonic.morphology_drift_radii[rate_index]
            perturbation = DiskMarginPerturbation(
                harmonic_index=harmonic.harmonic_index,
                rate_id=rate_id,
                center_shift_abs_upper=shift_float,
                old_response_radius=harmonic.response_radii[rate_index],
                new_response_radius=radius,
                old_morphology_radius=morphology_radius,
                new_morphology_radius=morphology_radius,
            )
            updates.append(update)
            perturbations.append(perturbation)
            update_proof.append(
                _update_payload(
                    update,
                    old_center,
                    harmonic.response_radii[rate_index],
                    morphology_radius,
                )
            )
            selected_manifest.append(
                {
                    "harmonic_index": harmonic.harmonic_index,
                    "rate_id": rate_id,
                    "window_id": link.request.window_id,
                    "window_result_sha256": canonical_window_result_sha256(row),
                    "record_manifest_id": row["record_manifest_id"],
                    "component_input_sha256": component_document["input_sha256"],
                    "component_certificate_sha256": component_document["certificate_sha256"],
                    "wls_execution_sha256": wls_evidence["execution_sha256"],
                    "wls_interval_certificate_sha256": wls_evidence[
                        "wls_interval_certificate_sha256"
                    ],
                    "external_digital_payload_sha256": wls_evidence[
                        "external_digital_payload_sha256"
                    ],
                    "target_exact_functional_radius_floor": wls_evidence[
                        "target_exact_functional_radius_floor"
                    ],
                }
            )
            replayed_wls_evidence.append(
                {
                    "harmonic_index": harmonic.harmonic_index,
                    "rate_id": rate_id,
                    **wls_evidence,
                }
            )
            embedded_components.append(
                {
                    "harmonic_index": harmonic.harmonic_index,
                    "rate_id": rate_id,
                    "certificate_json": component_json,
                }
            )

    selected_keys = {(value["window_id"], value["harmonic_index"]) for value in selected_manifest}
    if selected_keys != set(rows_by_key):
        return _fail(safe_manifest, "protocol_result_set_not_exactly_bound")

    response_disks = tuple(updates)
    try:
        updated_problem = _updated_problem(problem, response_disks)
        transfer = certify_margin_transfer(
            problem,
            box,
            baseline_certificate_json,
            tuple(perturbations),
        )
        transfer_replay = replay_margin_transfer_certificate(
            problem,
            box,
            baseline_certificate_json,
            tuple(perturbations),
            transfer.to_json(),
        )
    except (ArithmeticError, OverflowError, TypeError, ValueError):
        return _fail(safe_manifest, "margin_transfer_execution_failed")
    if not transfer_replay.valid:
        return _fail(safe_manifest, "margin_transfer_certificate_invalid")
    if transfer.relation not in {
        MarginTransferRelation.PRESERVED_ROBUST,
        MarginTransferRelation.PRESERVED_CLOSED,
    }:
        return _fail(safe_manifest, "margin_transfer_not_preserved")

    exact_sample_wls_covered = all(
        evidence["mathematical_exact_sample_wls_coverage"] is True
        for evidence in replayed_wls_evidence
    )
    if not exact_sample_wls_covered:
        relation = R2T1BundleRelation.NOT_CERTIFIABLE
        reason = "wls_exact_sample_accumulation_not_enclosed"
    elif transfer.relation is MarginTransferRelation.PRESERVED_ROBUST:
        relation = R2T1BundleRelation.PRESERVED_ROBUST
        reason = "exact_sample_wls_and_r2_family_preserve_robust_t1_witness"
    else:
        relation = R2T1BundleRelation.PRESERVED_CLOSED
        reason = "exact_sample_wls_and_r2_family_preserve_closed_t1_witness"

    family_sha = _sha256_json(update_proof)
    input_manifest = {
        "semantics": (
            "exact-sample-wls-centres-plus-r2-disks-preserve-frozen-old-t1-selector"
            if exact_sample_wls_covered
            else (
                "conditional-declared-r2-disks-with-replayed-binary64-wls-centres;"
                "high-level-preservation-requires-exact-sample-wls-enclosure"
            )
        ),
        "mathematical_exact_sample_wls_coverage": exact_sample_wls_covered,
        "uncovered_numeric_term": (
            None if exact_sample_wls_covered else "G_and_b_accumulation_rounding"
        ),
        "protocol_sha256": protocol_hash,
        "protocol_status": protocol_validation.status.value,
        "baseline_t1_certificate_json_sha256": _sha256_text(baseline_certificate_json),
        "baseline_t1_certificate_sha256": baseline_digest,
        "baseline_t1_input_sha256": baseline_input_hash,
        "selected_results_and_components": selected_manifest,
        "new_response_family_sha256": family_sha,
    }
    proof = {
        "protocol_validation": protocol_validation.as_dict(),
        "identity": {
            "subject_pseudonym": result_rows[0]["subject_pseudonym"],
            "lead_id": result_rows[0]["lead_id"],
            "physical_gain_id": result_rows[0]["physical_gain_id"],
            "signal_unit": _protocol_table(protocol_config, "units")["signal_unit"],
            "calibration_split_id": _protocol_table(protocol_config, "data_boundary")[
                "calibration_split_id"
            ],
            "identification_split_id": _protocol_table(protocol_config, "data_boundary")[
                "identification_split_id"
            ],
        },
        (
            "certified_response_disks"
            if exact_sample_wls_covered
            else "conditional_declared_response_disks"
        ): update_proof,
        "embedded_evidence": {
            "baseline_t1_certificate_json": baseline_certificate_json,
            "r2_component_certificates": embedded_components,
            "replayed_wls_execution_evidence": replayed_wls_evidence,
            "margin_transfer_certificate_json": transfer.to_json(),
        },
        (
            "certified_margin_transfer_relation"
            if exact_sample_wls_covered
            else "conditional_declared_disk_margin_transfer_relation"
        ): transfer.relation.value,
        "high_level_relation": relation.value,
        "mathematical_exact_sample_wls_coverage": exact_sample_wls_covered,
        "uncovered_numeric_term": (
            None if exact_sample_wls_covered else "G_and_b_accumulation_rounding"
        ),
    }
    return _Evaluation(
        relation,
        reason,
        input_manifest,
        proof,
        response_disks if exact_sample_wls_covered else (),
        updated_problem if exact_sample_wls_covered else None,
    )


def certify_r2_t1_bundle(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    protocol_config: Mapping[str, Any],
    baseline_certificate_json: str,
    links: Sequence[R2DiskCertificateLink],
) -> R2T1BundleCertificate:
    """Seal the full deterministic R2-to-T1 evidence chain."""

    fallback_manifest = _safe_manifest(
        problem,
        box,
        protocol_config,
        baseline_certificate_json,
        links,
    )
    try:
        evaluation = _evaluate(
            problem,
            box,
            protocol_config,
            baseline_certificate_json,
            links,
        )
    except (AttributeError, KeyError, OverflowError, TypeError, ValueError, RecursionError):
        evaluation = _fail(fallback_manifest, "bundle_input_invalid")
    evaluation = _apply_bundle_size_limit(evaluation, fallback_manifest)
    payload = _evaluation_payload(evaluation)
    certificate_json, digest = _seal(payload)
    return R2T1BundleCertificate(
        evaluation.relation,
        evaluation.reason,
        payload["input_sha256"],
        digest,
        certificate_json,
        evaluation.response_disks,
        evaluation.updated_problem,
    )


def replay_r2_t1_bundle(
    problem: CertifiedJointProblem,
    box: CertifiedParameterBox,
    protocol_config: Mapping[str, Any],
    baseline_certificate_json: str,
    links: Sequence[R2DiskCertificateLink],
    certificate: R2T1BundleCertificate | str,
) -> R2T1BundleReplayResult:
    """Recompute every child proof and reject any changed or missing link."""

    value = certificate.to_json() if isinstance(certificate, R2T1BundleCertificate) else certificate
    try:
        document = _load_json_object(value, byte_limit=MAX_CERTIFICATE_BYTES)
        digest = document.get("certificate_sha256")
        if not isinstance(digest, str) or len(digest) != 64:
            raise ValueError("invalid bundle digest")
        unsigned = dict(document)
        del unsigned["certificate_sha256"]
        if _sha256_json(unsigned) != digest:
            raise ValueError("bundle digest mismatch")
        if document.get("schema") != SCHEMA or document.get("algorithm") != ALGORITHM:
            raise ValueError("bundle schema mismatch")
        stored_relation = R2T1BundleRelation(document.get("relation"))
        if not isinstance(document.get("reason"), str) or not document["reason"]:
            raise ValueError("bundle reason invalid")
        if document.get("input_sha256") != _sha256_json(document.get("input_manifest")):
            raise ValueError("bundle input hash mismatch")
    except (KeyError, TypeError, ValueError, json.JSONDecodeError):
        return R2T1BundleReplayResult(
            False,
            R2T1BundleRelation.NOT_CERTIFIABLE,
            "invalid_bundle_envelope",
        )

    try:
        evaluation = _evaluate(
            problem,
            box,
            protocol_config,
            baseline_certificate_json,
            links,
        )
    except (AttributeError, KeyError, OverflowError, TypeError, ValueError, RecursionError):
        return R2T1BundleReplayResult(
            False,
            R2T1BundleRelation.NOT_CERTIFIABLE,
            "bundle_recomputation_failed",
        )
    fallback_manifest = _safe_manifest(
        problem,
        box,
        protocol_config,
        baseline_certificate_json,
        links,
    )
    evaluation = _apply_bundle_size_limit(evaluation, fallback_manifest)
    expected_proof = evaluation.proof
    if (
        stored_relation is not evaluation.relation
        or document.get("reason") != evaluation.reason
        or document.get("input_manifest") != evaluation.input_manifest
        or document.get("input_sha256") != _sha256_json(evaluation.input_manifest)
        or document.get("proof") != expected_proof
    ):
        return R2T1BundleReplayResult(
            False,
            R2T1BundleRelation.NOT_CERTIFIABLE,
            "bundle_proof_mismatch",
        )
    replay_reason = (
        "replay_verified"
        if evaluation.relation is not R2T1BundleRelation.NOT_CERTIFIABLE
        else "replay_not_certifiable"
    )
    return R2T1BundleReplayResult(
        True,
        evaluation.relation,
        replay_reason,
        evaluation.response_disks,
        evaluation.updated_problem,
    )


__all__ = [
    "ALGORITHM",
    "SCHEMA",
    "R2DiskCertificateLink",
    "R2T1BundleCertificate",
    "R2T1BundleRelation",
    "R2T1BundleReplayResult",
    "ResponseDiskUpdate",
    "certify_r2_t1_bundle",
    "replay_r2_t1_bundle",
]
