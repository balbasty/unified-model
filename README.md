# Unified model of Computational Anatomy

This is an attempt at unifying in a single Bayesian model different components classically used in computational anatomy.

* In its current (work in progress) form, this model performs:
  - Image segmentation
  - Missing data (either due to FOV or entirely missing modalities) imputation
  - Manual labels fusion
  - Image registration (i.e., normalisation to a standard space)
  - Template building
  - Shape analysis

* It deals with multichannel data (such as multi-contrast MRIs) and multiple modalities (e.g., CT and MRI).

* At the moment, all modalities of a single subject must have the same lattice (same voxel size, same dimensions). However, different subjects may have different lattices.

# Usage

The input of the algorithm is a `job` object, which can be a JSON file, with fields:
- `input.dat`: List of JSON files (or directories containing JSON files) that semantically describe all input volumes.
- [optional] `input.model`: List of JSON files (or directories containing JSON files) that describe some known parts of the model.
- `opt`: A (possibly hierarchical) structure of options.
- `par`: A structure of options related to parallel or distributed processing. This relies on the [`distributed-computing` toolbox](https://github.com/WTCN-computational-anatomy-group/distributed-computing).

## Input data

**At the moment, only unzipped Nifti files are handled by the software.**

All input data must be described in one or several JSON files. A typical JSON file will describe one or several input volumes.
If it describes one volume, it should contain a dictionary. If it describes several volumes, it should contain a list of dictionaries.

* Each dictionary should have the fields:
  - `name`: the subject name or ID.
  - [optional] `population`: the population (e.g., cohort) to which the subject belongs. If a population is provided, the subject unique identifier will be `<population>_<name>`.
  - [optional] `pth`: An absolute or relative (from the JSON file location) path to the described volume. If no path is provided, we will try to locate a Nifti file with the exact same name as the JSON file.

* If the described volume is a modality, the dictionary should have the field `modality` containing the modality name (e.g., `'CT'`).

* If the descibed volume is one channel from a multichannel modality, it should additionnaly have the field `channel` containing the channel name (e.g., `'T1'`).

* If the described volume is a manual segmentation, it should have the field  `rater` containing the rater ID and an optional field `protocol` containing the labelling protocol ID.

### Example

Here is an example file describing multichannel MRI data (2 subjects, 3 volumes / subject):

```python
[{'name':       'subj01',
  'population': 'MYDATABASE',
  'modality':   'MRI',
  'channel':    'T1,
  'pth':        'subj01_T1.nii',},
 {'name':       'subj01',
  'population': 'MYDATABASE',
  'modality':   'MRI',
  'channel':    'T2,
  'pth':        'subj01_T2.nii',},
 {'name':       'subj01',
  'population': 'MYDATABASE',
  'modality':   'MRI',
  'channel':    'PD,
  'pth':        'subj01_PD.nii',},
 {'name':       'subj02',
  'population': 'MYDATABASE',
  'modality':   'MRI',
  'channel':    'T1,
  'pth':        'subj02_T1.nii',},
 {'name':       'subj02',
  'population': 'MYDATABASE',
  'modality':   'MRI',
  'channel':    'T2,
  'pth':        'subj02_T2.nii',},
 {'name':       'subj02',
  'population': 'MYDATABASE',
  'modality':   'MRI',
  'channel':    'PD,
  'pth':        'subj02_PD.nii',},
]
```
Note that each dictionary could have been stored in a single JSON file (that could be named `'subj01_T1.json'`, etc.)

## Input model

**TODO**

## Options

One may want some options to be shared by the whole algorithm (e.g., such that all bias fields have the regularisation value `bias.prm = 0.1`) and some options to be specific to some parts only (e.g., only such that no bias field is estimated in CT images, `<if modality.name == 'CT'> bias.do = false `).

We tried to make this kind of modular option setting easy and user-friendly:
- to set a generic option, one would use (in JSON syntax):
  ```python
  'opt': {
    'bias': {'prm': 0.1}
  } 
  ```
- to set a specific option, one would use:
  ```python
  'opt': {
    'modality': [{
      'name': 'CT', 
      'bias': {'do': false}
    }]
  }
  ```
  
Here is the list of all possible options, in Matlab format. When an option can be set in both a generic and specific way, the optional hierarchical fields are written between squared brackets. Converting from Matlab to Python/JSON is quite straightforward (`a.b{i}.c` gives `['a']['b'](i)['c']`)
  
### Template

| Matlab                        | Default value   | Description | Choices |
|-------------------------------|-----------------|-------------|---------|
| `opt.template.do`   | `true`          | Optimise template?
| `opt.template.K`    | `6`             | Number of classes
| `opt.template.prm`  | `[1e-3 1e-1 0]` | Regularisation parameters | [Absolute Membrane Bending]
| `opt.template.bnd`  | `0`             | Boundary condition        | `0`=Circulant, `1`=Neumann
| `opt.template.itrp` | `1`             | Interpolation order
| `opt.template.vs`   | `1.5`           | Voxel size                | `NaN`=mean of input
| `opt.template.lat`  | `NaN`           | Lattice dimensions        | `NaN`=auto
| `opt.template.crop` | `true`          | Crop at BBox
| `opt.template.mrf`  | `false`         | Markov Random Field cleanup
| `opt.template.gn`   | `1`             | Number of Gauss-Newton updates

### Non-linear registration

| Matlab                         | Default value   | Description | Choices |
|--------------------------------|-----------------|-------------|---------|
| `opt.register.mode`  | `'diffeo'`      | Registration mode  | `'none'`, `'smalldef'`, `'diffeo'`, `'shape'`
| `opt.register.prm`   | `[0.001 0 10 0.1 0.2]` | Regularisation parameters | [Absolute Membrane Bending LE1 LE2]
| `opt.register.bnd`   | `0`             | Boundary condition | `0`=Circulant, `1`=Neumann
| `opt.register.gn`    | `1`             | Number of Gauss-Newton updates
| `opt.register.ls`    | `6`             | Number of Line-Search steps
| `opt.register.itg`   | `NaN`           | Number of shooting integration steps | `NaN`=auto
| `opt.register.write` | `true`          | Write deformation field (at the end)

### Shape

| Matlab                         | Default value   | Description | Choices |
|--------------------------------|-----------------|-------------|---------|
| `opt.register.shape.do` | `true`            | Optimise subpsace?
| `opt.register.shape.K`  | `64`              | Number of principal modes
| `opt.register.shape.A0` | `eye(K)`          | Prior latent precision matrix
| `opt.register.shape.n0` | `K`               | Prior latent degrees of freedom | `0`=ML, `inf`=Fixed
| `opt.register.shape.l0` | `17`              | Prior residual precision
| `opt.register.shape.m0` | `20`              | Prior residual degrees of freedom | `0`=ML, `inf`=Fixed
| `opt.register.shape.orthogonalise` | `true` | Orthogonalise subspace

### Affine

| Matlab                         | Default value   | Description | Choices |
|--------------------------------|-----------------|-------------|---------|
| `opt.register.affine.mode` | `'affine'` | Affine mode  | `'none'`, `'rigid'`, `' similitude'`, `'affine'`
| `opt.register.affine.gn`   | `1`        | Number of Gauss-Newton updates
| `opt.register.affine.ls`   | `6`        | Number of Line-Search steps
| `opt.register.affine.A0`   | `eye(Mr)`  | Prior precision matrix
| `opt.register.affine.n0`   | `Mr`       | Prior degrees of freedom

- **Note 1:** if `register.mode == 'shape'` then `affine.mode == 'rigid'`
- **Note 2:** `Mr` is the number of non-rigid parameters

### Bias correction

| Matlab                         | Default value   | Description | Choices |
|--------------------------------|-----------------|-------------|---------|
| `opt[.modality{m}][.channel{c}].bias.do`     | `true` | Perform bias correction?
| `opt[.modality{m}][.channel{c}].bias.prm`    | `1e-4` | Bending regularisation
| `opt[.modality{m}][.channel{c}].bias.gn`     | `1`    | Number of Gauss-Newton updates
| `opt[.modality{m}][.channel{c}].bias.ls`     | `6`    | Number of Line-Search steps
| `opt[.modality{m}][.channel{c}].bias.init`   | `'median'` | Method to initialse global scaling | `'median'`, `'mean'`
| `opt[.modality{m}][.channel{c}].bias.val0`   | `8192` | Value to match with the mean or median
| `opt[.modality{m}][.channel{c}].bias.center` | `true` | Zero-Center bias fields
| `opt[.modality{m}][.channel{c}].bias.write`  | `true` | Write bias field (at the end) 

- **Note**: if `modality.name == 'CT'` then `bias.do == false`

### Segmentation (Gaussian Mixture)

| Matlab                         | Default value   | Description | Choices |
|--------------------------------|-----------------|-------------|---------|
| `opt[.modality{m}].gmm.do`               | `true`       | Optimise intensity priors
| `opt[.modality{m}].gmm.K`                | `Ka`         | Number of Gaussian clusers
| `opt[.modality{m}].gmm.mix`              | `eye(Ka,K)`  | Map tissues to clusters (+ mixing dirichlet prior)
| `opt[.modality{m}].gmm.missing`          | `true`       | Infer missing values
| `opt[.modality{m}].gmm.constr`           | `false`      | Constrained precision matrices
| `opt[.modality{m}].gmm.iter`             | `20`         | Number of sub-EM iterations
| `opt[.modality{m}].gmm.init`             | `'kmeans'`   | Cluster initialisation | `'kmeans'`, `'gmm'`
| `opt[.modality{m}].gmm.order`            | `'magnitude'`| Cluster ordering | `'random'`, `'total'`, `'magnitude'`, `'user'`
| `opt[.modality{m}].gmm.V0`               | `eye(K)`     | Prior hyper-precision (if reg)
| `opt[.modality{m}].gmm.p0`               | `K`          | Prior hyper-precision degrees of freedom
| `opt[.modality{m}].gmm[.cluster{k}].mu0` | `'kmeans'`   | Prior mean
| `opt[.modality{m}].gmm[.cluster{k}].b0`  | `K`          | Prior mean degrees of freedom
| `opt[.modality{m}].gmm[.cluster{k}].A0`  | `'kmeans'`   | Prior precision
| `opt[.modality{m}].gmm[.cluster{k}].n0`  | `K`          | Prior precision degrees of freedom

- **Note:** `Ka` is the number of tissue classes in the template.

### General algorithm

| Matlab                         | Default value   | Description | Choices |
|--------------------------------|-----------------|-------------|---------|
| `opt.iter.max` | `30`   | Maximum number of VEM iterations
| `opt.iter.min` | `1`    | Minimum number of VEM iterations
| `opt.iter.tol` | `1e-4` | Lower-bound tolerance (stop if gain < tol)
| `opt.sample`   | `0`    | Sampling distance in mm to accelerate processing | `0`=no sampling

### User interface

| Matlab                         | Default value   | Description | Choices |
|--------------------------------|-----------------|-------------|---------|
| `opt.ui.verbose`      | `true`           | Print stuff during processing | `0`=nothing, `1`=normal, `2`=debug
| `opt.ui.figure.lb`    | `'lower bound'`  | Lower-bound tracking | `''`=do not plot
| `opt.ui.figure.model` | `'model'`        | Visualise model (template/shape/intensity) | `''`=do not plot
| `opt.ui.figure.reg`   | `'registration'` | Visualise registration (nat/warped/def/shape) | `''`=do not plot
| `opt.ui.figure.seg`   | `'segmentation'` | Visualise segmentation (int/nat/bias) | `''`=do not plot

## Future developments

In the future, this model could possibly include and support:
- Super-resolution / denoising
- Inputs with different lattices for a single subject
- Quantitative mapping
- ...

## Contributors

This software was developed by Mikael Brudfors and Yaël Balbastre in John Ashburner's [Computational Anatomy Group](http://www.fil.ion.ucl.ac.uk/Ashburner/) at the [Wellcome Centre for Human Neuroimaging](http://www.fil.ion.ucl.ac.uk/) in UCL.

If you encounter any difficulty, please send an email to `y.balbastre` or `mikael.brudfors.15` _at_ `ucl.ac.uk`

## License

This software is released under the [GNU General Public License version 3](LICENSE) (GPL v3). As a result, you may copy, distribute and modify the software as long as you track changes/dates in source files. Any modifications to or software including (via compiler) GPL-licensed code must also be made available under the GPL along with build & install instructions.

[TL;DR: GPL v3](https://tldrlegal.com/license/gnu-general-public-license-v3-(gpl-3))
