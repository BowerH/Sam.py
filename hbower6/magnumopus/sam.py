import re
class Read:
    def __init__(self, sam_line):
        fields = sam_line.strip().split('\t')
        self.qname = fields[0]
        self.flag = int(fields[1])
        self.rname = fields[2]
        self.pos = int(fields[3])
        self.mapq = int(fields[4])
        self.cigar = fields[5]
        self.rnext = fields[6]
        self.pnext = int(fields[7])
        self.tlen = int(fields[8])
        self.seq = fields[9]
        self.qual = fields[10]

    @property
    def is_mapped(self):
        return not (self.flag & 0x4)

    @property
    def is_forward(self):
        return not (self.flag & 0x10)

    @property
    def is_reverse(self):
        return bool(self.flag & 0x10)

    @property
    def is_primary(self):
        return not (self.flag & 0x100) and not (self.flag & 0x800)

    def mapped_seq(self) -> str:
        mapped = []
        read_position = 0
        for op, length in self._parse_cigar():
            if op in 'MX':
                start = read_position
                end = read_position + length
                segment = self.seq[start:end]
                mapped.append(segment)
                read_position += length
            elif op == 'I':
                start = read_position
                end = read_position + length
                inserted_segment = self.seq[start:end]
                mapped.append(inserted_segment)
                read_position += length
            elif op == 'D':
                deletion = '-' * length
                mapped.append(deletion)
            elif op == 'S':
                read_position += length
        return ''.join(mapped)

    def base_at_pos(self, ref_pos: int) -> str:
        if not self.is_mapped or ref_pos < self.pos:
            return ""

        current_ref_position = self.pos
        current_read_position = 0

        for op, (length) in self._parse_cigar():
            if op in ('MX'):
                if current_ref_position <= ref_pos < current_ref_position + length:
                    read_offset = ref_pos - current_ref_position
                    base_index = current_read_position + read_offset
                    base = self.seq[base_index]
                    return base
                current_ref_position += length
                current_read_position += length
            elif op == 'I':
                if current_ref_position == ref_pos:
                    start_index = current_read_position
                    end_index = current_read_position + length
                    inserted_bases = self.seq[start_index:end_index]
                    return inserted_bases
                current_read_position += length
            elif op == 'D':
                if current_ref_position <= ref_pos < current_ref_position + length:
                    return ""
                current_ref_position += length
            elif op == 'S':
                current_read_position += length

            if current_ref_position > ref_pos:
                break
        return ""

    def qual_at_pos(self, pos: int) -> str:
        if not self.is_mapped or pos < self.pos:
            return ""

        current_ref_position = self.pos
        current_read_position = 0

        for op, length in self._parse_cigar():
            if op in ('MX'):
                if current_ref_position <= pos < current_ref_position + length:
                    read_offset = pos - current_ref_position
                    quality_index = current_read_position + read_offset
                    quality_value = self.qual[quality_index]
                    return quality_value
                current_ref_position += length
                current_read_position += length
            elif op == 'I':
                current_read_position += length
            elif op == 'D':
                if current_ref_position <= pos < current_ref_position + length:
                    return ""
                current_ref_position += length
            elif op == 'S':
                current_read_position += length
        return ""

    def _parse_cigar(self):
        result = []
        length = ""
        for char in self.cigar:
            if char.isdigit():
                length += char
            else:
                result.append((char, int(length)))
                length = ""
        return result
