from collections import OrderedDict

def get_genome_sizes(genome_fai_path):
    lengths = {}
    for line in genome_fai_path.splitlines():
        if not line.strip():
            continue
        items = line.split()
        chrom = items[0]
        length = int(items[1])
        lengths[chrom] = length
    return lengths


class GenomeCoordinateConverter:
    def __init__(self, chrom_lens=None):
        if chrom_lens is None:
            chrom_lens = get_genome_sizes()

        sorted_chroms = sorted(chrom_lens.keys())
        offsets = OrderedDict()
        offset = 0
        chrom_spans = {}
        for chrom in sorted_chroms:
            chrom_len = chrom_lens[chrom]
            chrom_name = chrom.lower()
            offsets[chrom_name] = offset
            chrom_start = offset
            offset += chrom_len
            chrom_end = offset
            self.genome_size = offset
            chrom_spans[chrom_name] = (chrom_start, chrom_end)
        self._offsets = offsets
        self.chrom_spans = chrom_spans
        self.chrom_lens = chrom_lens

    def transform_coordinate(self, chrom, pos):
        if isinstance(chrom, bytes):
            chrom = chrom.decode()
        try:
            return self._offsets[chrom.lower()] + pos
        except KeyError:
            raise ValueError("Unknown chromosome: " + chrom)


class PositionInPericentromericRegion(Exception):
    pass


class GenomeCoordinateConverter2:
    def __init__(self):
        res = get_euchromatic_regions()
        self._euchromatic_regions = res["euchromatic_regions"]
        self.chrom_lens = res["euchromatic_sizes"]

        self.sorted_chroms = sorted(self._euchromatic_regions.keys())
        euchromatic_regions = self._euchromatic_regions

        chrom_offsets = {}
        chrom_spans = {}
        pericentromeric_starts = {}
        offset = 0
        for chrom in self.sorted_chroms:
            chrom_offsets[chrom] = offset

            chrom_start = offset
            offset += (
                euchromatic_regions[chrom][0][1]
                + euchromatic_regions[chrom][1][1]
                - euchromatic_regions[chrom][1][0]
            )
            pericentromeric_starts[chrom] = euchromatic_regions[chrom][0][1]
            chrom_end = offset
            self.genome_size = offset
            chrom_spans[chrom] = (chrom_start, chrom_end)

        self.chrom_offsets = chrom_offsets
        self.chrom_spans = chrom_spans
        self.pericentromeric_starts = pericentromeric_starts

    def _get_offset(self, chrom, pos):
        try:
            chrom_offset = self.chrom_offsets[chrom]
        except KeyError:
            raise KeyError("Unknown chromosome")

        euchromatic_regions = self._euchromatic_regions[chrom]
        if pos <= euchromatic_regions[0][1]:
            return chrom_offset, 0
        if euchromatic_regions[0][1] < pos and pos < euchromatic_regions[1][0]:
            raise PositionInPericentromericRegion()
        if pos >= euchromatic_regions[1][0]:
            return chrom_offset + euchromatic_regions[0][1], euchromatic_regions[1][0]

    def transform_coordinate(self, chrom, pos):
        if isinstance(chrom, bytes):
            chrom = chrom.decode()
        offset, to_remove_from_pos = self._get_offset(chrom, pos)
        return offset + pos - to_remove_from_pos