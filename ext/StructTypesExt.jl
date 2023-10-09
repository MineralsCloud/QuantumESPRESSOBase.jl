module StructTypesExt

using QuantumESPRESSOBase.PWscf:
    ControlNamelist,
    SystemNamelist,
    ElectronsNamelist,
    IonsNamelist,
    CellNamelist,
    DosNamelist,
    BandsNamelist,
    PWInput
using QuantumESPRESSOBase.PHonon: PhNamelist, Q2rNamelist, MatdynNamelist, DynmatNamelist
using StructTypes: Struct

import StructTypes: StructType

StructType(::Type{ControlNamelist}) = Struct()
StructType(::Type{SystemNamelist}) = Struct()
StructType(::Type{ElectronsNamelist}) = Struct()
StructType(::Type{IonsNamelist}) = Struct()
StructType(::Type{CellNamelist}) = Struct()
StructType(::Type{DosNamelist}) = Struct()
StructType(::Type{BandsNamelist}) = Struct()
StructType(::Type{PhNamelist}) = Struct()
StructType(::Type{Q2rNamelist}) = Struct()
StructType(::Type{MatdynNamelist}) = Struct()
StructType(::Type{DynmatNamelist}) = Struct()

end
