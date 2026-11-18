import gzip

import pytest

from vrs_anvil.coverage import vcf_to_bed

HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def row(pos, alt="T", chrom="chr1"):
    return f"{chrom}\t{pos}\t.\tA\t{alt}\t.\tLowQual\tEND=900\n"


@pytest.mark.parametrize("compressed", [False, True])
def test_union(tmp_path, compressed):
    source, output = tmp_path / "in.vcf", tmp_path / "out.bed"
    content = HEADER + row(100) + row(1) + row(101) + row(100) + row(5, chrom="chr2")
    source.write_bytes(
        gzip.compress(content.encode()) if compressed else content.encode()
    )
    vcf_to_bed(source, output, window=10)
    assert output.read_text() == "chr1\t0\t11\nchr1\t89\t111\nchr2\t0\t15\n"


@pytest.mark.parametrize("alt", ["T", "<DEL>", "A]chr2:100]", "<NON_REF>,T", "T,G"])
def test_variant_types(tmp_path, alt):
    source, output = tmp_path / "in", tmp_path / "out"
    source.write_text(HEADER + row(100, alt))
    vcf_to_bed(source, output)
    assert output.read_text() == "chr1\t99\t100\n"


@pytest.mark.parametrize("alt", [".", "*", "<NON_REF>", "<*>", "*,<NON_REF>"])
def test_placeholders(tmp_path, alt):
    source, output = tmp_path / "in", tmp_path / "out"
    source.write_text(HEADER + row(100, alt))
    vcf_to_bed(source, output)
    assert output.read_text() == ""


@pytest.mark.parametrize(
    "content", ["", HEADER + row(0), HEADER + row("bad"), HEADER + "chr1\t2\n"]
)
def test_invalid_preserves_output(tmp_path, content):
    source, output = tmp_path / "in", tmp_path / "out"
    source.write_text(content)
    output.write_text("existing")
    with pytest.raises(ValueError, match="Invalid VCF"):
        vcf_to_bed(source, output)
    assert output.read_text() == "existing"


@pytest.mark.parametrize("window", [-1, True, 1.5])
def test_invalid_window(tmp_path, window):
    with pytest.raises(ValueError, match="window"):
        vcf_to_bed(tmp_path / "in", tmp_path / "out", window=window)


def test_same_path(tmp_path):
    source = tmp_path / "in"
    source.write_text(HEADER)
    with pytest.raises(ValueError, match="paths must differ"):
        vcf_to_bed(source, source)
