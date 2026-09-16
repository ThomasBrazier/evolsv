"""A VCF writer that survives having no records to write.

Imported by vcf_sv_specification.py and light_vcf.py, both of which run under
the pysam pinned in workflow/envs/pysam_v2.yaml. Those scripts are invoked as
`python workflow/scripts/<name>.py`, so this module sits on sys.path[0] and
needs no packaging.
"""

from pysam import VariantFile


class VcfWriter:
    """Write records to a VCF, tolerating the case where there are none.

    pysam 0.10 (the version pinned in workflow/envs/pysam_v2.yaml) writes the
    VCF header lazily, on the first call to write(). Closing a writer that
    never received a record therefore dereferences an uninitialised header and
    segfaults the interpreter (exit status 139) instead of raising. That is how
    an input VCF holding a single variant used to kill vcf_sv_specification.py:
    its only variant was diverted to the `.ignored` file by the
    chromosome-edge filter, so the main writer was closed empty.

    The underlying pysam writer is therefore opened only once there is
    something to write, and the header is emitted directly when there is not --
    the downstream `bcftools view` still needs a well-formed, if variant-free,
    VCF rather than the zero-byte file pysam leaves behind.
    """

    def __init__(self, path, header):
        self.path = path
        self.header = header
        self._out = None
        self._closed = False

    def write(self, record):
        if self._out is None:
            self._out = VariantFile(self.path, "w", header=self.header)
        self._out.write(record)

    def close(self):
        if self._closed:
            return
        self._closed = True
        if self._out is None:
            with open(self.path, "w") as empty:
                empty.write(str(self.header))
        else:
            self._out.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()
