using DataStructures

using DataStructures
helsinki_eeg_bad_channels = DefaultDict((key) -> (@warn "Bad channels undefined for patient $key"; []), Dict(
    9 => ["ECG","Resp","Cz"],
    19 => ["ECG","Resp","Cz","Fz"], #unverified
    21 => ["ECG","Resp","Cz","Fz"], #unverified
    31 => ["ECG","Resp","Cz"],
    44 => ["ECG","Resp","Cz"],
    47 => ["ECG","Resp","Cz","Fz"],
    50 => ["ECG","Resp","Cz"],
    62 => ["ECG","Resp","Cz","Fz"], #unverified
    75 => ["ECG","Resp","Cz","Fz"]
); passkey=true)