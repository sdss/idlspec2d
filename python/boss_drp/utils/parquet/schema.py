from pydl.pydlutils.yanny import read_table_yanny
import pyarrow as pa


# ----------------------------
# Build Arrow schema from Yanny datamodel
# ----------------------------

class Schema:
    def __init__(self, yanny_file=None):
        self.yanny_file = yanny_file
        self.schema_def = None
        self.column_meta = {}
        self.primary_hdr = {}
        self.mapping = None
        self.numeric_cols = []
        self.id_cols = []

    def datamodel_arrow_schema(self, struct_name="EXT1"):
        rows = read_table_yanny(self.yanny_file, struct_name)
        rows.convert_bytestring_to_unicode()

        fields = []
        for row in rows:
            name = row["Column"]
            ytype = row["type"]
            unit = row['unit'] 
            cnull = row['null']
            #TODO: add ucd 
            
            description = row['description']
            if ytype[:-1].isdigit():
                n = int(ytype[:-1])
                base = ytype[-1]
            else:
                n = None
                base = ytype

            if base == "K":
                arrow_type = pa.int64()
            elif base == "I" or base == "J":
                arrow_type = pa.int32()
            elif base == "B":
                arrow_type = pa.int64()
            elif base == "E":
                arrow_type = pa.float32()
            elif base == "D":
                arrow_type = pa.float64()
            elif base == "A":
                arrow_type = pa.string()
            else:
                raise ValueError(f"Unknown type {ytype} for {name}")

            if name == 'SDSS5_TARGET_FLAGS':
                arrow_type = pa.list_(arrow_type)
            elif n is not None:
                arrow_type = pa.list_(arrow_type, n)


            fields.append(pa.field(name, arrow_type))

            self.column_meta[name] = {
                "unit": unit,
                "description": description,
                "fits_type": ytype,
                "null": cnull
            }
        self.schema_def = pa.schema(fields)

    def datamodel_header_metadata(self, struct_name="HDR0"):
        rows = read_table_yanny(self.yanny_file, struct_name)
        rows.convert_bytestring_to_unicode()

        for row in rows:
            key = row['card']
            description = row["description"]
            self.primary_hdr[key] = description

    def get_mapping(self, struct_name="MAPPING"):
        rows = read_table_yanny(self.yanny_file, struct_name)
        rows.convert_bytestring_to_unicode()

        self.numeric_cols  = rows[rows['numeric'] == 1]['spAllColumn'].tolist()
        self.id_cols = rows[rows['numeric'] == 0]['spAllColumn'].tolist()
        self.mapping = rows

