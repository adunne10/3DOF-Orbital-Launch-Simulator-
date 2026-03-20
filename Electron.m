
function vehicle = Electron()

vehicle.name = "Electron";

 vehicle.stage  = [
 struct( ...
    'm_dry',   572, ...
    'm_prop',  9314, ...
    'thrust',  229164, ...
    'isp',     311, ...
    'Cd',      0.2, ...
    'A',       1.13, ...
    'A_exit',  0.0491, ...
    'pld_guess', 100 ) %lower bounds

struct( ...
    'm_dry',   254, ...
    'm_prop',  2860, ...
    'thrust',  25800, ...
    'isp',     343, ...
    'Cd',      0.2, ...
    'A',       1.13, ...
    'A_exit',  0 , ...
    'pld_guess', 400) %upper bounds
];
vehicle.performance.interval=100e3;

vehicle.performance.altitude=[
    400000
    500000
    600000
    700000
    800000
    900000
    1000000
    1100000
    1200000
    ];

vehicle.performance.payload_40=[
    270
    265
    255
    250
    245
    237
    230
    222
    215
    ];

vehicle.performance.payload_60=[
    250
    242
    238
    231
    222
    218
    211
    204
    198
    ];

vehicle.performance.payload_80=[
    222
    219
    212
    208
    201
    195
    189
    185
    176
    ];

vehicle.performance.payload_SSO=[
    202
    198
    192
    187
    180
    175
    168
    160
    156
    ];

end

