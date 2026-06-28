# nano-wgs rerun instructions

This workflow is run with:

```bash
bash nano-wgs.snakejob_rerun_28062026.sh \
    metadata_patients_rerun_28062026.txt \
    /data/cephfs-2/unmirrored/projects/henssen-nanopore \
    Evolution
```

The arguments are:

```text
metadata_file = metadata_patients_rerun_28062026.txt
workdir       = /data/cephfs-2/unmirrored/projects/henssen-nanopore
project       = Evolution
```

## Input path format

The workflow expects raw FASTQ files to be located under:

```text
{workdir}/{project}/{run}/{subdir}
```

In the Snakefile, the FASTQ input directories are constructed as:

```python
fastq=lambda wildcards: [
    "{workdir}/{project}/{run}/{subdir}".format(
        workdir=wildcards.workdir,
        project=wildcards.project,
        run=row.Run,
        subdir=NAME_FASTQ_FOLDER
    )
    for index, row in metadata[metadata.Sample == wildcards.sample].iterrows()
]
```

Therefore, the `Run` column in the metadata file must contain the path **relative to `{workdir}/{project}`**, without the final `pass` folder.

In this case:

```python
NAME_FASTQ_FOLDER = "pass"
```

So the full input path is built as:

```text
{workdir}/{project}/{Run}/pass
```

## Example with `fastq_guppy634`

If the full FASTQ directory is:

```text
/data/cephfs-2/unmirrored/projects/henssen-nanopore/Evolution/RawData/CB1008-2nd-look_01012020/CB1008-2nd-look/CB1008-2nd-look-1/20210421_1245_MN31855_FAP90140_c703d4a8/fastq_guppy634/pass
```

then the metadata `Run` value should be:

```text
RawData/CB1008-2nd-look_01012020/CB1008-2nd-look/CB1008-2nd-look-1/20210421_1245_MN31855_FAP90140_c703d4a8/fastq_guppy634
```

The workflow then adds:

```text
pass
```

automatically.

## Example with `fastq`

Another valid metadata `Run` value is:

```text
RawData/CB1008-2nd-look_01012020/CB1008-2nd-look/CB1008-2nd-look-1/20210421_1245_MN31855_FAP90140_c703d4a8/fastq
```

In that case, the workflow will search in:

```text
/data/cephfs-2/unmirrored/projects/henssen-nanopore/Evolution/RawData/CB1008-2nd-look_01012020/CB1008-2nd-look/CB1008-2nd-look-1/20210421_1245_MN31855_FAP90140_c703d4a8/fastq/pass
```

## Important check before running

Before running the workflow, make sure that the constructed input path exists.

For example:

```bash
ls -lh /data/cephfs-2/unmirrored/projects/henssen-nanopore/Evolution/RawData/CB1008-2nd-look_01012020/CB1008-2nd-look/CB1008-2nd-look-1/20210421_1245_MN31855_FAP90140_c703d4a8/fastq_guppy634/pass
```

The key point is:

```text
Full input path = workdir + project + Run + NAME_FASTQ_FOLDER
```

For this rerun:

```text
/data/cephfs-2/unmirrored/projects/henssen-nanopore
+ Evolution
+ RawData/.../fastq_guppy634
+ pass
```

which becomes:

```text
/data/cephfs-2/unmirrored/projects/henssen-nanopore/Evolution/RawData/.../fastq_guppy634/pass
```

If Snakemake raises a `MissingInputException`, first check that the `Run` value in the metadata file is correct and that the final constructed path exists on disk.
