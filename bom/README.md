# BOM — Compressor Test Rig

Bill of Materials, organised by component type. Each subsystem domain has its
own YAML file. `bom_index.yaml` is the master registry used by aggregation tooling.

---

## File Structure

```
bom/
├── bom_index.yaml                     ← master registry (start here)
│
├── bom_mechanical_rotating.yaml       ← rotor, shaft, bearings, seals, tie-bolts
├── bom_mechanical_static.yaml         ← bellmouth, ogive centerbody, casing, diffuser
├── bom_structure.yaml                 ← bedplate, pedestals, motor base, guarding
├── bom_fasteners.yaml                 ← general structural bolts, dowels
│
├── bom_drivetrain.yaml                ← main motor, VFD, flexible coupling
│
├── bom_instrumentation_pressure.yaml  ← Scanivalve scanner, taps, transducers
├── bom_instrumentation_flow.yaml      ← bellmouth dP, traverse probes, hot-wire
├── bom_instrumentation_temp.yaml      ← thermocouples, RTDs, signal conditioners
├── bom_instrumentation_vibration.yaml ← accelerometers, proximity probes, DAQ
│
├── bom_motion_control.yaml            ← HIWIN stage, AM8012, AX5101, ballscrew
│
├── bom_electrical_power.yaml          ← MCC, breakers, EStop, power cables
├── bom_electrical_signal.yaml         ← EtherCAT, thermocouple cables, Ethernet
├── bom_controls_hardware.yaml         ← Beckhoff IPC, EtherCAT terminals, I/O
│
├── bom_pneumatics_hydraulics.yaml     ← lube oil system, seal air (TBD)
└── bom_consumables.yaml               ← gaskets, O-rings, lubricants, sealants
```

---

## Entry Schema

Every component entry uses these fields:

| Field | Description |
|---|---|
| `part_number` | Internal PN (`CTR-XXX-NNN`) or manufacturer PN |
| `description` | Plain-language description |
| `qty` | Quantity per rig build |
| `unit` | `ea` / `m` / `kg` / `set` |
| `manufacturer` | OEM name |
| `supplier` | Distributor / purchasing contact |
| `unit_cost_eur` | Estimated unit cost in EUR |
| `lead_time_weeks` | Procurement lead time |
| `drawing_number` | Associated design drawing reference |
| `material` | Material specification |
| `notes` | Free-text remarks, cross-references, warnings |
| `status` | `tbd` \| `quoted` \| `ordered` \| `received` \| `installed` |

Rotating parts additionally carry a `balance_record` field pointing to the
dynamic balance report in `docs/balance_records/`.

---

## Status Workflow

```
tbd → quoted → ordered → received → installed
```

Update `status` in each entry as procurement progresses. Update `revision` in
`bom_index.yaml` on any structural change to a sub-BOM file.

---

## Scope Boundaries

A few non-obvious decisions on what goes where:

- **Rotor tie-bolts** → `bom_mechanical_rotating.yaml` (torque-critical rotor hardware, not general fasteners)
- **Instrument port fittings / Swagelok** → respective instrumentation file
- **AX5101 servo drive** → `bom_controls_hardware.yaml`; **AM8012 motor** → `bom_motion_control.yaml`
- **Main drive VFD** → `bom_drivetrain.yaml` (not controls hardware)
- **EtherCAT cabling** → `bom_electrical_signal.yaml`; **EtherCAT terminals** → `bom_controls_hardware.yaml`

---

## Tooling

`tools/aggregate_bom.py` reads `bom_index.yaml`, loads all sub-BOM files, and
produces `bom_master.csv` and `bom_master.xlsx` with a cost roll-up and
procurement status summary. *(Script TBD.)*

---

## Key Cross-References

| If you are specifying... | Also check... |
|---|---|
| Bearings (rotating) | Pedestal bores in `bom_structure.yaml` |
| Motor (drivetrain) | Motor support foot pattern in `bom_structure.yaml` |
| Thermocouples | Beckhoff EL3314/EL3204 channel count in `bom_controls_hardware.yaml` |
| Scanivalve scanner | Ethernet port on IPC in `bom_controls_hardware.yaml` |
| Proximity probes | Probe boss provisions in `bom_structure.yaml` |
| Throttle valve (actuated) | Analogue output terminals in `bom_controls_hardware.yaml` |
| Traverse probes | Linear stage in `bom_motion_control.yaml` |
