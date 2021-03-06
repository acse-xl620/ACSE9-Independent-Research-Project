??
??
D
AddV2
x"T
y"T
z"T"
Ttype:
2	??
B
AssignVariableOp
resource
value"dtype"
dtypetype?
~
BiasAdd

value"T	
bias"T
output"T" 
Ttype:
2	"-
data_formatstringNHWC:
NHWCNCHW
8
Const
output"dtype"
valuetensor"
dtypetype
?
Conv2D

input"T
filter"T
output"T"
Ttype:	
2"
strides	list(int)"
use_cudnn_on_gpubool(",
paddingstring:
SAMEVALIDEXPLICIT""
explicit_paddings	list(int)
 "-
data_formatstringNHWC:
NHWCNCHW" 
	dilations	list(int)

;
Elu
features"T
activations"T"
Ttype:
2
?
FusedBatchNormV3
x"T

scale"U
offset"U	
mean"U
variance"U
y"T

batch_mean"U
batch_variance"U
reserve_space_1"U
reserve_space_2"U
reserve_space_3"U"
Ttype:
2"
Utype:
2"
epsilonfloat%??8"&
exponential_avg_factorfloat%  ??";
data_formatstringNHWC:
NHWCNCHWNDHWCNCDHW"
is_trainingbool(
.
Identity

input"T
output"T"	
Ttype
q
MatMul
a"T
b"T
product"T"
transpose_abool( "
transpose_bbool( "
Ttype:

2	
e
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool(?
?
Mul
x"T
y"T
z"T"
Ttype:
2	?

NoOp
M
Pack
values"T*N
output"T"
Nint(0"	
Ttype"
axisint 
C
Placeholder
output"dtype"
dtypetype"
shapeshape:
@
ReadVariableOp
resource
value"dtype"
dtypetype?
[
Reshape
tensor"T
shape"Tshape
output"T"	
Ttype"
Tshapetype0:
2	
?
ResizeNearestNeighbor
images"T
size
resized_images"T"
Ttype:
2
	"
align_cornersbool( "
half_pixel_centersbool( 
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0?
.
Rsqrt
x"T
y"T"
Ttype:

2
l
SaveV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0?
?
Select
	condition

t"T
e"T
output"T"	
Ttype
P
Shape

input"T
output"out_type"	
Ttype"
out_typetype0:
2	
H
ShardedFilename
basename	
shard

num_shards
filename
0
Sigmoid
x"T
y"T"
Ttype:

2
?
StatefulPartitionedCall
args2Tin
output2Tout"
Tin
list(type)("
Tout
list(type)("	
ffunc"
configstring "
config_protostring "
executor_typestring ?
@
StaticRegexFullMatch	
input

output
"
patternstring
?
StridedSlice

input"T
begin"Index
end"Index
strides"Index
output"T"	
Ttype"
Indextype:
2	"

begin_maskint "
end_maskint "
ellipsis_maskint "
new_axis_maskint "
shrink_axis_maskint 
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
<
Sub
x"T
y"T
z"T"
Ttype:
2	
?
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 ?"serve*2.6.02v2.6.0-0-g919f693420e8??
y
dense_2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	d?*
shared_namedense_2/kernel
r
"dense_2/kernel/Read/ReadVariableOpReadVariableOpdense_2/kernel*
_output_shapes
:	d?*
dtype0
q
dense_2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:?*
shared_namedense_2/bias
j
 dense_2/bias/Read/ReadVariableOpReadVariableOpdense_2/bias*
_output_shapes	
:?*
dtype0
?
batch_normalization_3/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:?*,
shared_namebatch_normalization_3/gamma
?
/batch_normalization_3/gamma/Read/ReadVariableOpReadVariableOpbatch_normalization_3/gamma*
_output_shapes	
:?*
dtype0
?
batch_normalization_3/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:?*+
shared_namebatch_normalization_3/beta
?
.batch_normalization_3/beta/Read/ReadVariableOpReadVariableOpbatch_normalization_3/beta*
_output_shapes	
:?*
dtype0
?
!batch_normalization_3/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:?*2
shared_name#!batch_normalization_3/moving_mean
?
5batch_normalization_3/moving_mean/Read/ReadVariableOpReadVariableOp!batch_normalization_3/moving_mean*
_output_shapes	
:?*
dtype0
?
%batch_normalization_3/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:?*6
shared_name'%batch_normalization_3/moving_variance
?
9batch_normalization_3/moving_variance/Read/ReadVariableOpReadVariableOp%batch_normalization_3/moving_variance*
_output_shapes	
:?*
dtype0
?
conv2d_2/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:@@* 
shared_nameconv2d_2/kernel
{
#conv2d_2/kernel/Read/ReadVariableOpReadVariableOpconv2d_2/kernel*&
_output_shapes
:@@*
dtype0
r
conv2d_2/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*
shared_nameconv2d_2/bias
k
!conv2d_2/bias/Read/ReadVariableOpReadVariableOpconv2d_2/bias*
_output_shapes
:@*
dtype0
?
batch_normalization_4/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*,
shared_namebatch_normalization_4/gamma
?
/batch_normalization_4/gamma/Read/ReadVariableOpReadVariableOpbatch_normalization_4/gamma*
_output_shapes
:@*
dtype0
?
batch_normalization_4/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*+
shared_namebatch_normalization_4/beta
?
.batch_normalization_4/beta/Read/ReadVariableOpReadVariableOpbatch_normalization_4/beta*
_output_shapes
:@*
dtype0
?
!batch_normalization_4/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*2
shared_name#!batch_normalization_4/moving_mean
?
5batch_normalization_4/moving_mean/Read/ReadVariableOpReadVariableOp!batch_normalization_4/moving_mean*
_output_shapes
:@*
dtype0
?
%batch_normalization_4/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:@*6
shared_name'%batch_normalization_4/moving_variance
?
9batch_normalization_4/moving_variance/Read/ReadVariableOpReadVariableOp%batch_normalization_4/moving_variance*
_output_shapes
:@*
dtype0
?
conv2d_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:@?* 
shared_nameconv2d_3/kernel
|
#conv2d_3/kernel/Read/ReadVariableOpReadVariableOpconv2d_3/kernel*'
_output_shapes
:@?*
dtype0
s
conv2d_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:?*
shared_nameconv2d_3/bias
l
!conv2d_3/bias/Read/ReadVariableOpReadVariableOpconv2d_3/bias*
_output_shapes	
:?*
dtype0
?
batch_normalization_5/gammaVarHandleOp*
_output_shapes
: *
dtype0*
shape:?*,
shared_namebatch_normalization_5/gamma
?
/batch_normalization_5/gamma/Read/ReadVariableOpReadVariableOpbatch_normalization_5/gamma*
_output_shapes	
:?*
dtype0
?
batch_normalization_5/betaVarHandleOp*
_output_shapes
: *
dtype0*
shape:?*+
shared_namebatch_normalization_5/beta
?
.batch_normalization_5/beta/Read/ReadVariableOpReadVariableOpbatch_normalization_5/beta*
_output_shapes	
:?*
dtype0
?
!batch_normalization_5/moving_meanVarHandleOp*
_output_shapes
: *
dtype0*
shape:?*2
shared_name#!batch_normalization_5/moving_mean
?
5batch_normalization_5/moving_mean/Read/ReadVariableOpReadVariableOp!batch_normalization_5/moving_mean*
_output_shapes	
:?*
dtype0
?
%batch_normalization_5/moving_varianceVarHandleOp*
_output_shapes
: *
dtype0*
shape:?*6
shared_name'%batch_normalization_5/moving_variance
?
9batch_normalization_5/moving_variance/Read/ReadVariableOpReadVariableOp%batch_normalization_5/moving_variance*
_output_shapes	
:?*
dtype0
y
dense_3/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	?@2*
shared_namedense_3/kernel
r
"dense_3/kernel/Read/ReadVariableOpReadVariableOpdense_3/kernel*
_output_shapes
:	?@2*
dtype0
p
dense_3/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:2*
shared_namedense_3/bias
i
 dense_3/bias/Read/ReadVariableOpReadVariableOpdense_3/bias*
_output_shapes
:2*
dtype0

NoOpNoOp
?7
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*?6
value?6B?6 B?6
?
layer_with_weights-0
layer-0
layer-1
layer_with_weights-1
layer-2
layer-3
layer_with_weights-2
layer-4
layer-5
layer_with_weights-3
layer-6
layer_with_weights-4
layer-7
	layer-8

layer_with_weights-5

layer-9
layer-10
layer_with_weights-6
layer-11
trainable_variables
	variables
regularization_losses
	keras_api

signatures
h

kernel
bias
trainable_variables
	variables
regularization_losses
	keras_api
R
trainable_variables
	variables
regularization_losses
	keras_api
?
axis
	gamma
beta
moving_mean
 moving_variance
!trainable_variables
"	variables
#regularization_losses
$	keras_api
R
%trainable_variables
&	variables
'regularization_losses
(	keras_api
h

)kernel
*bias
+trainable_variables
,	variables
-regularization_losses
.	keras_api
R
/trainable_variables
0	variables
1regularization_losses
2	keras_api
?
3axis
	4gamma
5beta
6moving_mean
7moving_variance
8trainable_variables
9	variables
:regularization_losses
;	keras_api
h

<kernel
=bias
>trainable_variables
?	variables
@regularization_losses
A	keras_api
R
Btrainable_variables
C	variables
Dregularization_losses
E	keras_api
?
Faxis
	Ggamma
Hbeta
Imoving_mean
Jmoving_variance
Ktrainable_variables
L	variables
Mregularization_losses
N	keras_api
R
Otrainable_variables
P	variables
Qregularization_losses
R	keras_api
h

Skernel
Tbias
Utrainable_variables
V	variables
Wregularization_losses
X	keras_api
f
0
1
2
3
)4
*5
46
57
<8
=9
G10
H11
S12
T13
?
0
1
2
3
4
 5
)6
*7
48
59
610
711
<12
=13
G14
H15
I16
J17
S18
T19
 
?
Ymetrics
trainable_variables
	variables

Zlayers
[layer_metrics
\layer_regularization_losses
regularization_losses
]non_trainable_variables
 
ZX
VARIABLE_VALUEdense_2/kernel6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_2/bias4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE

0
1

0
1
 
?
^metrics
trainable_variables
_non_trainable_variables
	variables
`layer_metrics
alayer_regularization_losses
regularization_losses

blayers
 
 
 
?
cmetrics
trainable_variables
dnon_trainable_variables
	variables
elayer_metrics
flayer_regularization_losses
regularization_losses

glayers
 
fd
VARIABLE_VALUEbatch_normalization_3/gamma5layer_with_weights-1/gamma/.ATTRIBUTES/VARIABLE_VALUE
db
VARIABLE_VALUEbatch_normalization_3/beta4layer_with_weights-1/beta/.ATTRIBUTES/VARIABLE_VALUE
rp
VARIABLE_VALUE!batch_normalization_3/moving_mean;layer_with_weights-1/moving_mean/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUE%batch_normalization_3/moving_variance?layer_with_weights-1/moving_variance/.ATTRIBUTES/VARIABLE_VALUE

0
1

0
1
2
 3
 
?
hmetrics
!trainable_variables
inon_trainable_variables
"	variables
jlayer_metrics
klayer_regularization_losses
#regularization_losses

llayers
 
 
 
?
mmetrics
%trainable_variables
nnon_trainable_variables
&	variables
olayer_metrics
player_regularization_losses
'regularization_losses

qlayers
[Y
VARIABLE_VALUEconv2d_2/kernel6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE
WU
VARIABLE_VALUEconv2d_2/bias4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE

)0
*1

)0
*1
 
?
rmetrics
+trainable_variables
snon_trainable_variables
,	variables
tlayer_metrics
ulayer_regularization_losses
-regularization_losses

vlayers
 
 
 
?
wmetrics
/trainable_variables
xnon_trainable_variables
0	variables
ylayer_metrics
zlayer_regularization_losses
1regularization_losses

{layers
 
fd
VARIABLE_VALUEbatch_normalization_4/gamma5layer_with_weights-3/gamma/.ATTRIBUTES/VARIABLE_VALUE
db
VARIABLE_VALUEbatch_normalization_4/beta4layer_with_weights-3/beta/.ATTRIBUTES/VARIABLE_VALUE
rp
VARIABLE_VALUE!batch_normalization_4/moving_mean;layer_with_weights-3/moving_mean/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUE%batch_normalization_4/moving_variance?layer_with_weights-3/moving_variance/.ATTRIBUTES/VARIABLE_VALUE

40
51

40
51
62
73
 
?
|metrics
8trainable_variables
}non_trainable_variables
9	variables
~layer_metrics
layer_regularization_losses
:regularization_losses
?layers
[Y
VARIABLE_VALUEconv2d_3/kernel6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUE
WU
VARIABLE_VALUEconv2d_3/bias4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUE

<0
=1

<0
=1
 
?
?metrics
>trainable_variables
?non_trainable_variables
?	variables
?layer_metrics
 ?layer_regularization_losses
@regularization_losses
?layers
 
 
 
?
?metrics
Btrainable_variables
?non_trainable_variables
C	variables
?layer_metrics
 ?layer_regularization_losses
Dregularization_losses
?layers
 
fd
VARIABLE_VALUEbatch_normalization_5/gamma5layer_with_weights-5/gamma/.ATTRIBUTES/VARIABLE_VALUE
db
VARIABLE_VALUEbatch_normalization_5/beta4layer_with_weights-5/beta/.ATTRIBUTES/VARIABLE_VALUE
rp
VARIABLE_VALUE!batch_normalization_5/moving_mean;layer_with_weights-5/moving_mean/.ATTRIBUTES/VARIABLE_VALUE
zx
VARIABLE_VALUE%batch_normalization_5/moving_variance?layer_with_weights-5/moving_variance/.ATTRIBUTES/VARIABLE_VALUE

G0
H1

G0
H1
I2
J3
 
?
?metrics
Ktrainable_variables
?non_trainable_variables
L	variables
?layer_metrics
 ?layer_regularization_losses
Mregularization_losses
?layers
 
 
 
?
?metrics
Otrainable_variables
?non_trainable_variables
P	variables
?layer_metrics
 ?layer_regularization_losses
Qregularization_losses
?layers
ZX
VARIABLE_VALUEdense_3/kernel6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_3/bias4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUE

S0
T1

S0
T1
 
?
?metrics
Utrainable_variables
?non_trainable_variables
V	variables
?layer_metrics
 ?layer_regularization_losses
Wregularization_losses
?layers
 
V
0
1
2
3
4
5
6
7
	8

9
10
11
 
 
*
0
 1
62
73
I4
J5
 
 
 
 
 
 
 
 
 
 
 

0
 1
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

60
71
 
 
 
 
 
 
 
 
 
 
 
 
 
 

I0
J1
 
 
 
 
 
 
 
 
 
 
 
 
 
?
serving_default_dense_2_inputPlaceholder*'
_output_shapes
:?????????d*
dtype0*
shape:?????????d
?
StatefulPartitionedCallStatefulPartitionedCallserving_default_dense_2_inputdense_2/kerneldense_2/bias%batch_normalization_3/moving_variancebatch_normalization_3/gamma!batch_normalization_3/moving_meanbatch_normalization_3/betaconv2d_2/kernelconv2d_2/biasbatch_normalization_4/gammabatch_normalization_4/beta!batch_normalization_4/moving_mean%batch_normalization_4/moving_varianceconv2d_3/kernelconv2d_3/biasbatch_normalization_5/gammabatch_normalization_5/beta!batch_normalization_5/moving_mean%batch_normalization_5/moving_variancedense_3/kerneldense_3/bias* 
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????2*6
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8? *.
f)R'
%__inference_signature_wrapper_5330683
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
?	
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename"dense_2/kernel/Read/ReadVariableOp dense_2/bias/Read/ReadVariableOp/batch_normalization_3/gamma/Read/ReadVariableOp.batch_normalization_3/beta/Read/ReadVariableOp5batch_normalization_3/moving_mean/Read/ReadVariableOp9batch_normalization_3/moving_variance/Read/ReadVariableOp#conv2d_2/kernel/Read/ReadVariableOp!conv2d_2/bias/Read/ReadVariableOp/batch_normalization_4/gamma/Read/ReadVariableOp.batch_normalization_4/beta/Read/ReadVariableOp5batch_normalization_4/moving_mean/Read/ReadVariableOp9batch_normalization_4/moving_variance/Read/ReadVariableOp#conv2d_3/kernel/Read/ReadVariableOp!conv2d_3/bias/Read/ReadVariableOp/batch_normalization_5/gamma/Read/ReadVariableOp.batch_normalization_5/beta/Read/ReadVariableOp5batch_normalization_5/moving_mean/Read/ReadVariableOp9batch_normalization_5/moving_variance/Read/ReadVariableOp"dense_3/kernel/Read/ReadVariableOp dense_3/bias/Read/ReadVariableOpConst*!
Tin
2*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *)
f$R"
 __inference__traced_save_5331876
?
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenamedense_2/kerneldense_2/biasbatch_normalization_3/gammabatch_normalization_3/beta!batch_normalization_3/moving_mean%batch_normalization_3/moving_varianceconv2d_2/kernelconv2d_2/biasbatch_normalization_4/gammabatch_normalization_4/beta!batch_normalization_4/moving_mean%batch_normalization_4/moving_varianceconv2d_3/kernelconv2d_3/biasbatch_normalization_5/gammabatch_normalization_5/beta!batch_normalization_5/moving_mean%batch_normalization_5/moving_variancedense_3/kerneldense_3/bias* 
Tin
2*
Tout
2*
_collective_manager_ids
 *
_output_shapes
: * 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *,
f'R%
#__inference__traced_restore_5331946??
?
?
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5331544

inputs&
readvariableop_resource:	?(
readvariableop_1_resource:	?7
(fusedbatchnormv3_readvariableop_resource:	?9
*fusedbatchnormv3_readvariableop_1_resource:	?
identity??FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1u
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp{
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes	
:?*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes	
:?*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*b
_output_shapesP
N:,????????????????????????????:?:?:?:?:*
epsilon%o?:*
is_training( 2
FusedBatchNormV3?
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*B
_output_shapes0
.:,????????????????????????????2

Identity?
NoOpNoOp ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:,????????????????????????????: : : : 2B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:j f
B
_output_shapes0
.:,????????????????????????????
 
_user_specified_nameinputs
?	
?
7__inference_batch_normalization_5_layer_call_fn_5331624

inputs
unknown:	?
	unknown_0:	?
	unknown_1:	?
	unknown_2:	?
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *B
_output_shapes0
.:,????????????????????????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_53295812
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*B
_output_shapes0
.:,????????????????????????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:,????????????????????????????: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:j f
B
_output_shapes0
.:,????????????????????????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5330017

inputs&
readvariableop_resource:	?(
readvariableop_1_resource:	?7
(fusedbatchnormv3_readvariableop_resource:	?9
*fusedbatchnormv3_readvariableop_1_resource:	?
identity??AssignNewValue?AssignNewValue_1?FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1u
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp{
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes	
:?*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes	
:?*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*P
_output_shapes>
<:??????????:?:?:?:?:*
epsilon%o?:*
exponential_avg_factor%
?#<2
FusedBatchNormV3?
AssignNewValueAssignVariableOp(fusedbatchnormv3_readvariableop_resourceFusedBatchNormV3:batch_mean:0 ^FusedBatchNormV3/ReadVariableOp*
_output_shapes
 *
dtype02
AssignNewValue?
AssignNewValue_1AssignVariableOp*fusedbatchnormv3_readvariableop_1_resource!FusedBatchNormV3:batch_variance:0"^FusedBatchNormV3/ReadVariableOp_1*
_output_shapes
 *
dtype02
AssignNewValue_1x
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*0
_output_shapes
:??????????2

Identity?
NoOpNoOp^AssignNewValue^AssignNewValue_1 ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:??????????: : : : 2 
AssignNewValueAssignNewValue2$
AssignNewValue_1AssignNewValue_12B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:X T
0
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
%__inference_signature_wrapper_5330683
dense_2_input
unknown:	d?
	unknown_0:	?
	unknown_1:	?
	unknown_2:	?
	unknown_3:	?
	unknown_4:	?#
	unknown_5:@@
	unknown_6:@
	unknown_7:@
	unknown_8:@
	unknown_9:@

unknown_10:@%

unknown_11:@?

unknown_12:	?

unknown_13:	?

unknown_14:	?

unknown_15:	?

unknown_16:	?

unknown_17:	?@2

unknown_18:2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCalldense_2_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18* 
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????2*6
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8? *+
f&R$
"__inference__wrapped_model_53291552
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????22

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*N
_input_shapes=
;:?????????d: : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:V R
'
_output_shapes
:?????????d
'
_user_specified_namedense_2_input
?
?
*__inference_conv2d_2_layer_call_fn_5331298

inputs!
unknown:@@
	unknown_0:@
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv2d_2_layer_call_and_return_conditional_losses_53297322
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:?????????@2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????@: : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
??
?
I__inference_sequential_1_layer_call_and_return_conditional_losses_5330828

inputs9
&dense_2_matmul_readvariableop_resource:	d?6
'dense_2_biasadd_readvariableop_resource:	?F
7batch_normalization_3_batchnorm_readvariableop_resource:	?J
;batch_normalization_3_batchnorm_mul_readvariableop_resource:	?H
9batch_normalization_3_batchnorm_readvariableop_1_resource:	?H
9batch_normalization_3_batchnorm_readvariableop_2_resource:	?A
'conv2d_2_conv2d_readvariableop_resource:@@6
(conv2d_2_biasadd_readvariableop_resource:@;
-batch_normalization_4_readvariableop_resource:@=
/batch_normalization_4_readvariableop_1_resource:@L
>batch_normalization_4_fusedbatchnormv3_readvariableop_resource:@N
@batch_normalization_4_fusedbatchnormv3_readvariableop_1_resource:@B
'conv2d_3_conv2d_readvariableop_resource:@?7
(conv2d_3_biasadd_readvariableop_resource:	?<
-batch_normalization_5_readvariableop_resource:	?>
/batch_normalization_5_readvariableop_1_resource:	?M
>batch_normalization_5_fusedbatchnormv3_readvariableop_resource:	?O
@batch_normalization_5_fusedbatchnormv3_readvariableop_1_resource:	?9
&dense_3_matmul_readvariableop_resource:	?@25
'dense_3_biasadd_readvariableop_resource:2
identity??.batch_normalization_3/batchnorm/ReadVariableOp?0batch_normalization_3/batchnorm/ReadVariableOp_1?0batch_normalization_3/batchnorm/ReadVariableOp_2?2batch_normalization_3/batchnorm/mul/ReadVariableOp?5batch_normalization_4/FusedBatchNormV3/ReadVariableOp?7batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1?$batch_normalization_4/ReadVariableOp?&batch_normalization_4/ReadVariableOp_1?5batch_normalization_5/FusedBatchNormV3/ReadVariableOp?7batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1?$batch_normalization_5/ReadVariableOp?&batch_normalization_5/ReadVariableOp_1?conv2d_2/BiasAdd/ReadVariableOp?conv2d_2/Conv2D/ReadVariableOp?/conv2d_2/bias/Regularizer/Square/ReadVariableOp?1conv2d_2/kernel/Regularizer/Square/ReadVariableOp?conv2d_3/BiasAdd/ReadVariableOp?conv2d_3/Conv2D/ReadVariableOp?/conv2d_3/bias/Regularizer/Square/ReadVariableOp?1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?dense_2/BiasAdd/ReadVariableOp?dense_2/MatMul/ReadVariableOp?.dense_2/bias/Regularizer/Square/ReadVariableOp?0dense_2/kernel/Regularizer/Square/ReadVariableOp?dense_3/BiasAdd/ReadVariableOp?dense_3/MatMul/ReadVariableOp?.dense_3/bias/Regularizer/Square/ReadVariableOp?0dense_3/kernel/Regularizer/Square/ReadVariableOp?
dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes
:	d?*
dtype02
dense_2/MatMul/ReadVariableOp?
dense_2/MatMulMatMulinputs%dense_2/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
dense_2/MatMul?
dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02 
dense_2/BiasAdd/ReadVariableOp?
dense_2/BiasAddBiasAdddense_2/MatMul:product:0&dense_2/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
dense_2/BiasAddn
dense_2/EluEludense_2/BiasAdd:output:0*
T0*(
_output_shapes
:??????????2
dense_2/Elu?
dropout_1/IdentityIdentitydense_2/Elu:activations:0*
T0*(
_output_shapes
:??????????2
dropout_1/Identity?
.batch_normalization_3/batchnorm/ReadVariableOpReadVariableOp7batch_normalization_3_batchnorm_readvariableop_resource*
_output_shapes	
:?*
dtype020
.batch_normalization_3/batchnorm/ReadVariableOp?
%batch_normalization_3/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2'
%batch_normalization_3/batchnorm/add/y?
#batch_normalization_3/batchnorm/addAddV26batch_normalization_3/batchnorm/ReadVariableOp:value:0.batch_normalization_3/batchnorm/add/y:output:0*
T0*
_output_shapes	
:?2%
#batch_normalization_3/batchnorm/add?
%batch_normalization_3/batchnorm/RsqrtRsqrt'batch_normalization_3/batchnorm/add:z:0*
T0*
_output_shapes	
:?2'
%batch_normalization_3/batchnorm/Rsqrt?
2batch_normalization_3/batchnorm/mul/ReadVariableOpReadVariableOp;batch_normalization_3_batchnorm_mul_readvariableop_resource*
_output_shapes	
:?*
dtype024
2batch_normalization_3/batchnorm/mul/ReadVariableOp?
#batch_normalization_3/batchnorm/mulMul)batch_normalization_3/batchnorm/Rsqrt:y:0:batch_normalization_3/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2%
#batch_normalization_3/batchnorm/mul?
%batch_normalization_3/batchnorm/mul_1Muldropout_1/Identity:output:0'batch_normalization_3/batchnorm/mul:z:0*
T0*(
_output_shapes
:??????????2'
%batch_normalization_3/batchnorm/mul_1?
0batch_normalization_3/batchnorm/ReadVariableOp_1ReadVariableOp9batch_normalization_3_batchnorm_readvariableop_1_resource*
_output_shapes	
:?*
dtype022
0batch_normalization_3/batchnorm/ReadVariableOp_1?
%batch_normalization_3/batchnorm/mul_2Mul8batch_normalization_3/batchnorm/ReadVariableOp_1:value:0'batch_normalization_3/batchnorm/mul:z:0*
T0*
_output_shapes	
:?2'
%batch_normalization_3/batchnorm/mul_2?
0batch_normalization_3/batchnorm/ReadVariableOp_2ReadVariableOp9batch_normalization_3_batchnorm_readvariableop_2_resource*
_output_shapes	
:?*
dtype022
0batch_normalization_3/batchnorm/ReadVariableOp_2?
#batch_normalization_3/batchnorm/subSub8batch_normalization_3/batchnorm/ReadVariableOp_2:value:0)batch_normalization_3/batchnorm/mul_2:z:0*
T0*
_output_shapes	
:?2%
#batch_normalization_3/batchnorm/sub?
%batch_normalization_3/batchnorm/add_1AddV2)batch_normalization_3/batchnorm/mul_1:z:0'batch_normalization_3/batchnorm/sub:z:0*
T0*(
_output_shapes
:??????????2'
%batch_normalization_3/batchnorm/add_1{
reshape_1/ShapeShape)batch_normalization_3/batchnorm/add_1:z:0*
T0*
_output_shapes
:2
reshape_1/Shape?
reshape_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
reshape_1/strided_slice/stack?
reshape_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2!
reshape_1/strided_slice/stack_1?
reshape_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2!
reshape_1/strided_slice/stack_2?
reshape_1/strided_sliceStridedSlicereshape_1/Shape:output:0&reshape_1/strided_slice/stack:output:0(reshape_1/strided_slice/stack_1:output:0(reshape_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
reshape_1/strided_slicex
reshape_1/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2
reshape_1/Reshape/shape/1x
reshape_1/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :2
reshape_1/Reshape/shape/2x
reshape_1/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :@2
reshape_1/Reshape/shape/3?
reshape_1/Reshape/shapePack reshape_1/strided_slice:output:0"reshape_1/Reshape/shape/1:output:0"reshape_1/Reshape/shape/2:output:0"reshape_1/Reshape/shape/3:output:0*
N*
T0*
_output_shapes
:2
reshape_1/Reshape/shape?
reshape_1/ReshapeReshape)batch_normalization_3/batchnorm/add_1:z:0 reshape_1/Reshape/shape:output:0*
T0*/
_output_shapes
:?????????@2
reshape_1/Reshape?
conv2d_2/Conv2D/ReadVariableOpReadVariableOp'conv2d_2_conv2d_readvariableop_resource*&
_output_shapes
:@@*
dtype02 
conv2d_2/Conv2D/ReadVariableOp?
conv2d_2/Conv2DConv2Dreshape_1/Reshape:output:0&conv2d_2/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????@*
paddingSAME*
strides
2
conv2d_2/Conv2D?
conv2d_2/BiasAdd/ReadVariableOpReadVariableOp(conv2d_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02!
conv2d_2/BiasAdd/ReadVariableOp?
conv2d_2/BiasAddBiasAddconv2d_2/Conv2D:output:0'conv2d_2/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????@2
conv2d_2/BiasAddx
conv2d_2/EluEluconv2d_2/BiasAdd:output:0*
T0*/
_output_shapes
:?????????@2
conv2d_2/Elu{
up_sampling2d/ConstConst*
_output_shapes
:*
dtype0*
valueB"      2
up_sampling2d/Const
up_sampling2d/Const_1Const*
_output_shapes
:*
dtype0*
valueB"      2
up_sampling2d/Const_1?
up_sampling2d/mulMulup_sampling2d/Const:output:0up_sampling2d/Const_1:output:0*
T0*
_output_shapes
:2
up_sampling2d/mul?
*up_sampling2d/resize/ResizeNearestNeighborResizeNearestNeighborconv2d_2/Elu:activations:0up_sampling2d/mul:z:0*
T0*/
_output_shapes
:?????????@*
half_pixel_centers(2,
*up_sampling2d/resize/ResizeNearestNeighbor?
$batch_normalization_4/ReadVariableOpReadVariableOp-batch_normalization_4_readvariableop_resource*
_output_shapes
:@*
dtype02&
$batch_normalization_4/ReadVariableOp?
&batch_normalization_4/ReadVariableOp_1ReadVariableOp/batch_normalization_4_readvariableop_1_resource*
_output_shapes
:@*
dtype02(
&batch_normalization_4/ReadVariableOp_1?
5batch_normalization_4/FusedBatchNormV3/ReadVariableOpReadVariableOp>batch_normalization_4_fusedbatchnormv3_readvariableop_resource*
_output_shapes
:@*
dtype027
5batch_normalization_4/FusedBatchNormV3/ReadVariableOp?
7batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1ReadVariableOp@batch_normalization_4_fusedbatchnormv3_readvariableop_1_resource*
_output_shapes
:@*
dtype029
7batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1?
&batch_normalization_4/FusedBatchNormV3FusedBatchNormV3;up_sampling2d/resize/ResizeNearestNeighbor:resized_images:0,batch_normalization_4/ReadVariableOp:value:0.batch_normalization_4/ReadVariableOp_1:value:0=batch_normalization_4/FusedBatchNormV3/ReadVariableOp:value:0?batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*K
_output_shapes9
7:?????????@:@:@:@:@:*
epsilon%o?:*
is_training( 2(
&batch_normalization_4/FusedBatchNormV3?
conv2d_3/Conv2D/ReadVariableOpReadVariableOp'conv2d_3_conv2d_readvariableop_resource*'
_output_shapes
:@?*
dtype02 
conv2d_3/Conv2D/ReadVariableOp?
conv2d_3/Conv2DConv2D*batch_normalization_4/FusedBatchNormV3:y:0&conv2d_3/Conv2D/ReadVariableOp:value:0*
T0*0
_output_shapes
:??????????*
paddingSAME*
strides
2
conv2d_3/Conv2D?
conv2d_3/BiasAdd/ReadVariableOpReadVariableOp(conv2d_3_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02!
conv2d_3/BiasAdd/ReadVariableOp?
conv2d_3/BiasAddBiasAddconv2d_3/Conv2D:output:0'conv2d_3/BiasAdd/ReadVariableOp:value:0*
T0*0
_output_shapes
:??????????2
conv2d_3/BiasAddy
conv2d_3/EluEluconv2d_3/BiasAdd:output:0*
T0*0
_output_shapes
:??????????2
conv2d_3/Elu
up_sampling2d_1/ConstConst*
_output_shapes
:*
dtype0*
valueB"      2
up_sampling2d_1/Const?
up_sampling2d_1/Const_1Const*
_output_shapes
:*
dtype0*
valueB"      2
up_sampling2d_1/Const_1?
up_sampling2d_1/mulMulup_sampling2d_1/Const:output:0 up_sampling2d_1/Const_1:output:0*
T0*
_output_shapes
:2
up_sampling2d_1/mul?
,up_sampling2d_1/resize/ResizeNearestNeighborResizeNearestNeighborconv2d_3/Elu:activations:0up_sampling2d_1/mul:z:0*
T0*0
_output_shapes
:??????????*
half_pixel_centers(2.
,up_sampling2d_1/resize/ResizeNearestNeighbor?
$batch_normalization_5/ReadVariableOpReadVariableOp-batch_normalization_5_readvariableop_resource*
_output_shapes	
:?*
dtype02&
$batch_normalization_5/ReadVariableOp?
&batch_normalization_5/ReadVariableOp_1ReadVariableOp/batch_normalization_5_readvariableop_1_resource*
_output_shapes	
:?*
dtype02(
&batch_normalization_5/ReadVariableOp_1?
5batch_normalization_5/FusedBatchNormV3/ReadVariableOpReadVariableOp>batch_normalization_5_fusedbatchnormv3_readvariableop_resource*
_output_shapes	
:?*
dtype027
5batch_normalization_5/FusedBatchNormV3/ReadVariableOp?
7batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1ReadVariableOp@batch_normalization_5_fusedbatchnormv3_readvariableop_1_resource*
_output_shapes	
:?*
dtype029
7batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1?
&batch_normalization_5/FusedBatchNormV3FusedBatchNormV3=up_sampling2d_1/resize/ResizeNearestNeighbor:resized_images:0,batch_normalization_5/ReadVariableOp:value:0.batch_normalization_5/ReadVariableOp_1:value:0=batch_normalization_5/FusedBatchNormV3/ReadVariableOp:value:0?batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*P
_output_shapes>
<:??????????:?:?:?:?:*
epsilon%o?:*
is_training( 2(
&batch_normalization_5/FusedBatchNormV3s
flatten_1/ConstConst*
_output_shapes
:*
dtype0*
valueB"????    2
flatten_1/Const?
flatten_1/ReshapeReshape*batch_normalization_5/FusedBatchNormV3:y:0flatten_1/Const:output:0*
T0*(
_output_shapes
:??????????@2
flatten_1/Reshape?
dense_3/MatMul/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes
:	?@2*
dtype02
dense_3/MatMul/ReadVariableOp?
dense_3/MatMulMatMulflatten_1/Reshape:output:0%dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????22
dense_3/MatMul?
dense_3/BiasAdd/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype02 
dense_3/BiasAdd/ReadVariableOp?
dense_3/BiasAddBiasAdddense_3/MatMul:product:0&dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????22
dense_3/BiasAddy
dense_3/SigmoidSigmoiddense_3/BiasAdd:output:0*
T0*'
_output_shapes
:?????????22
dense_3/Sigmoid?
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes
:	d?*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp?
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	d?2#
!dense_2/kernel/Regularizer/Square?
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const?
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/Sum?
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_2/kernel/Regularizer/mul/x?
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mul?
.dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype020
.dense_2/bias/Regularizer/Square/ReadVariableOp?
dense_2/bias/Regularizer/SquareSquare6dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2!
dense_2/bias/Regularizer/Square?
dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_2/bias/Regularizer/Const?
dense_2/bias/Regularizer/SumSum#dense_2/bias/Regularizer/Square:y:0'dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/Sum?
dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_2/bias/Regularizer/mul/x?
dense_2/bias/Regularizer/mulMul'dense_2/bias/Regularizer/mul/x:output:0%dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/mul?
1conv2d_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'conv2d_2_conv2d_readvariableop_resource*&
_output_shapes
:@@*
dtype023
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_2/kernel/Regularizer/SquareSquare9conv2d_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*&
_output_shapes
:@@2$
"conv2d_2/kernel/Regularizer/Square?
!conv2d_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_2/kernel/Regularizer/Const?
conv2d_2/kernel/Regularizer/SumSum&conv2d_2/kernel/Regularizer/Square:y:0*conv2d_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/Sum?
!conv2d_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_2/kernel/Regularizer/mul/x?
conv2d_2/kernel/Regularizer/mulMul*conv2d_2/kernel/Regularizer/mul/x:output:0(conv2d_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/mul?
/conv2d_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp(conv2d_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype021
/conv2d_2/bias/Regularizer/Square/ReadVariableOp?
 conv2d_2/bias/Regularizer/SquareSquare7conv2d_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2"
 conv2d_2/bias/Regularizer/Square?
conv2d_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_2/bias/Regularizer/Const?
conv2d_2/bias/Regularizer/SumSum$conv2d_2/bias/Regularizer/Square:y:0(conv2d_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/Sum?
conv2d_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_2/bias/Regularizer/mul/x?
conv2d_2/bias/Regularizer/mulMul(conv2d_2/bias/Regularizer/mul/x:output:0&conv2d_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/mul?
1conv2d_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'conv2d_3_conv2d_readvariableop_resource*'
_output_shapes
:@?*
dtype023
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_3/kernel/Regularizer/SquareSquare9conv2d_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*'
_output_shapes
:@?2$
"conv2d_3/kernel/Regularizer/Square?
!conv2d_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_3/kernel/Regularizer/Const?
conv2d_3/kernel/Regularizer/SumSum&conv2d_3/kernel/Regularizer/Square:y:0*conv2d_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/Sum?
!conv2d_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_3/kernel/Regularizer/mul/x?
conv2d_3/kernel/Regularizer/mulMul*conv2d_3/kernel/Regularizer/mul/x:output:0(conv2d_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/mul?
/conv2d_3/bias/Regularizer/Square/ReadVariableOpReadVariableOp(conv2d_3_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype021
/conv2d_3/bias/Regularizer/Square/ReadVariableOp?
 conv2d_3/bias/Regularizer/SquareSquare7conv2d_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2"
 conv2d_3/bias/Regularizer/Square?
conv2d_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_3/bias/Regularizer/Const?
conv2d_3/bias/Regularizer/SumSum$conv2d_3/bias/Regularizer/Square:y:0(conv2d_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/Sum?
conv2d_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_3/bias/Regularizer/mul/x?
conv2d_3/bias/Regularizer/mulMul(conv2d_3/bias/Regularizer/mul/x:output:0&conv2d_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/mul?
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes
:	?@2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp?
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?@22#
!dense_3/kernel/Regularizer/Square?
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const?
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum?
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_3/kernel/Regularizer/mul/x?
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul?
.dense_3/bias/Regularizer/Square/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype020
.dense_3/bias/Regularizer/Square/ReadVariableOp?
dense_3/bias/Regularizer/SquareSquare6dense_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:22!
dense_3/bias/Regularizer/Square?
dense_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_3/bias/Regularizer/Const?
dense_3/bias/Regularizer/SumSum#dense_3/bias/Regularizer/Square:y:0'dense_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/Sum?
dense_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_3/bias/Regularizer/mul/x?
dense_3/bias/Regularizer/mulMul'dense_3/bias/Regularizer/mul/x:output:0%dense_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/muln
IdentityIdentitydense_3/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:?????????22

Identity?

NoOpNoOp/^batch_normalization_3/batchnorm/ReadVariableOp1^batch_normalization_3/batchnorm/ReadVariableOp_11^batch_normalization_3/batchnorm/ReadVariableOp_23^batch_normalization_3/batchnorm/mul/ReadVariableOp6^batch_normalization_4/FusedBatchNormV3/ReadVariableOp8^batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1%^batch_normalization_4/ReadVariableOp'^batch_normalization_4/ReadVariableOp_16^batch_normalization_5/FusedBatchNormV3/ReadVariableOp8^batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1%^batch_normalization_5/ReadVariableOp'^batch_normalization_5/ReadVariableOp_1 ^conv2d_2/BiasAdd/ReadVariableOp^conv2d_2/Conv2D/ReadVariableOp0^conv2d_2/bias/Regularizer/Square/ReadVariableOp2^conv2d_2/kernel/Regularizer/Square/ReadVariableOp ^conv2d_3/BiasAdd/ReadVariableOp^conv2d_3/Conv2D/ReadVariableOp0^conv2d_3/bias/Regularizer/Square/ReadVariableOp2^conv2d_3/kernel/Regularizer/Square/ReadVariableOp^dense_2/BiasAdd/ReadVariableOp^dense_2/MatMul/ReadVariableOp/^dense_2/bias/Regularizer/Square/ReadVariableOp1^dense_2/kernel/Regularizer/Square/ReadVariableOp^dense_3/BiasAdd/ReadVariableOp^dense_3/MatMul/ReadVariableOp/^dense_3/bias/Regularizer/Square/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*N
_input_shapes=
;:?????????d: : : : : : : : : : : : : : : : : : : : 2`
.batch_normalization_3/batchnorm/ReadVariableOp.batch_normalization_3/batchnorm/ReadVariableOp2d
0batch_normalization_3/batchnorm/ReadVariableOp_10batch_normalization_3/batchnorm/ReadVariableOp_12d
0batch_normalization_3/batchnorm/ReadVariableOp_20batch_normalization_3/batchnorm/ReadVariableOp_22h
2batch_normalization_3/batchnorm/mul/ReadVariableOp2batch_normalization_3/batchnorm/mul/ReadVariableOp2n
5batch_normalization_4/FusedBatchNormV3/ReadVariableOp5batch_normalization_4/FusedBatchNormV3/ReadVariableOp2r
7batch_normalization_4/FusedBatchNormV3/ReadVariableOp_17batch_normalization_4/FusedBatchNormV3/ReadVariableOp_12L
$batch_normalization_4/ReadVariableOp$batch_normalization_4/ReadVariableOp2P
&batch_normalization_4/ReadVariableOp_1&batch_normalization_4/ReadVariableOp_12n
5batch_normalization_5/FusedBatchNormV3/ReadVariableOp5batch_normalization_5/FusedBatchNormV3/ReadVariableOp2r
7batch_normalization_5/FusedBatchNormV3/ReadVariableOp_17batch_normalization_5/FusedBatchNormV3/ReadVariableOp_12L
$batch_normalization_5/ReadVariableOp$batch_normalization_5/ReadVariableOp2P
&batch_normalization_5/ReadVariableOp_1&batch_normalization_5/ReadVariableOp_12B
conv2d_2/BiasAdd/ReadVariableOpconv2d_2/BiasAdd/ReadVariableOp2@
conv2d_2/Conv2D/ReadVariableOpconv2d_2/Conv2D/ReadVariableOp2b
/conv2d_2/bias/Regularizer/Square/ReadVariableOp/conv2d_2/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp1conv2d_2/kernel/Regularizer/Square/ReadVariableOp2B
conv2d_3/BiasAdd/ReadVariableOpconv2d_3/BiasAdd/ReadVariableOp2@
conv2d_3/Conv2D/ReadVariableOpconv2d_3/Conv2D/ReadVariableOp2b
/conv2d_3/bias/Regularizer/Square/ReadVariableOp/conv2d_3/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp1conv2d_3/kernel/Regularizer/Square/ReadVariableOp2@
dense_2/BiasAdd/ReadVariableOpdense_2/BiasAdd/ReadVariableOp2>
dense_2/MatMul/ReadVariableOpdense_2/MatMul/ReadVariableOp2`
.dense_2/bias/Regularizer/Square/ReadVariableOp.dense_2/bias/Regularizer/Square/ReadVariableOp2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp2@
dense_3/BiasAdd/ReadVariableOpdense_3/BiasAdd/ReadVariableOp2>
dense_3/MatMul/ReadVariableOpdense_3/MatMul/ReadVariableOp2`
.dense_3/bias/Regularizer/Square/ReadVariableOp.dense_3/bias/Regularizer/Square/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
e
F__inference_dropout_1_layer_call_and_return_conditional_losses_5331145

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Constt
dropout/MulMulinputsdropout/Const:output:0*
T0*(
_output_shapes
:??????????2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*(
_output_shapes
:??????????*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:??????????2
dropout/GreaterEqual?
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:??????????2
dropout/Cast{
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*(
_output_shapes
:??????????2
dropout/Mul_1f
IdentityIdentitydropout/Mul_1:z:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:??????????:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
D__inference_dense_2_layer_call_and_return_conditional_losses_5329671

inputs1
matmul_readvariableop_resource:	d?.
biasadd_readvariableop_resource:	?
identity??BiasAdd/ReadVariableOp?MatMul/ReadVariableOp?.dense_2/bias/Regularizer/Square/ReadVariableOp?0dense_2/kernel/Regularizer/Square/ReadVariableOp?
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	d?*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2	
BiasAddV
EluEluBiasAdd:output:0*
T0*(
_output_shapes
:??????????2
Elu?
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	d?*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp?
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	d?2#
!dense_2/kernel/Regularizer/Square?
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const?
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/Sum?
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_2/kernel/Regularizer/mul/x?
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mul?
.dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:?*
dtype020
.dense_2/bias/Regularizer/Square/ReadVariableOp?
dense_2/bias/Regularizer/SquareSquare6dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2!
dense_2/bias/Regularizer/Square?
dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_2/bias/Regularizer/Const?
dense_2/bias/Regularizer/SumSum#dense_2/bias/Regularizer/Square:y:0'dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/Sum?
dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_2/bias/Regularizer/mul/x?
dense_2/bias/Regularizer/mulMul'dense_2/bias/Regularizer/mul/x:output:0%dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/mulm
IdentityIdentityElu:activations:0^NoOp*
T0*(
_output_shapes
:??????????2

Identity?
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp/^dense_2/bias/Regularizer/Square/ReadVariableOp1^dense_2/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:?????????d: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2`
.dense_2/bias/Regularizer/Square/ReadVariableOp.dense_2/bias/Regularizer/Square/ReadVariableOp2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
d
+__inference_dropout_1_layer_call_fn_5331155

inputs
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_dropout_1_layer_call_and_return_conditional_losses_53301282
StatefulPartitionedCall|
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:??????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:??????????22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5330076

inputs%
readvariableop_resource:@'
readvariableop_1_resource:@6
(fusedbatchnormv3_readvariableop_resource:@8
*fusedbatchnormv3_readvariableop_1_resource:@
identity??AssignNewValue?AssignNewValue_1?FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1t
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
:@*
dtype02
ReadVariableOpz
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes
:@*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes
:@*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes
:@*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*K
_output_shapes9
7:?????????@:@:@:@:@:*
epsilon%o?:*
exponential_avg_factor%
?#<2
FusedBatchNormV3?
AssignNewValueAssignVariableOp(fusedbatchnormv3_readvariableop_resourceFusedBatchNormV3:batch_mean:0 ^FusedBatchNormV3/ReadVariableOp*
_output_shapes
 *
dtype02
AssignNewValue?
AssignNewValue_1AssignVariableOp*fusedbatchnormv3_readvariableop_1_resource!FusedBatchNormV3:batch_variance:0"^FusedBatchNormV3/ReadVariableOp_1*
_output_shapes
 *
dtype02
AssignNewValue_1w
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*/
_output_shapes
:?????????@2

Identity?
NoOpNoOp^AssignNewValue^AssignNewValue_1 ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*6
_input_shapes%
#:?????????@: : : : 2 
AssignNewValueAssignNewValue2$
AssignNewValue_1AssignNewValue_12B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5331364

inputs%
readvariableop_resource:@'
readvariableop_1_resource:@6
(fusedbatchnormv3_readvariableop_resource:@8
*fusedbatchnormv3_readvariableop_1_resource:@
identity??AssignNewValue?AssignNewValue_1?FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1t
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
:@*
dtype02
ReadVariableOpz
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes
:@*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes
:@*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes
:@*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*]
_output_shapesK
I:+???????????????????????????@:@:@:@:@:*
epsilon%o?:*
exponential_avg_factor%
?#<2
FusedBatchNormV3?
AssignNewValueAssignVariableOp(fusedbatchnormv3_readvariableop_resourceFusedBatchNormV3:batch_mean:0 ^FusedBatchNormV3/ReadVariableOp*
_output_shapes
 *
dtype02
AssignNewValue?
AssignNewValue_1AssignVariableOp*fusedbatchnormv3_readvariableop_1_resource!FusedBatchNormV3:batch_variance:0"^FusedBatchNormV3/ReadVariableOp_1*
_output_shapes
 *
dtype02
AssignNewValue_1?
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*A
_output_shapes/
-:+???????????????????????????@2

Identity?
NoOpNoOp^AssignNewValue^AssignNewValue_1 ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*H
_input_shapes7
5:+???????????????????????????@: : : : 2 
AssignNewValueAssignNewValue2$
AssignNewValue_1AssignNewValue_12B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:i e
A
_output_shapes/
-:+???????????????????????????@
 
_user_specified_nameinputs
?
M
1__inference_up_sampling2d_1_layer_call_fn_5331526

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_up_sampling2d_1_layer_call_and_return_conditional_losses_53298102
PartitionedCallu
IdentityIdentityPartitionedCall:output:0*
T0*0
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:??????????:X T
0
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
)__inference_dense_2_layer_call_fn_5331128

inputs
unknown:	d?
	unknown_0:	?
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_2_layer_call_and_return_conditional_losses_53296712
StatefulPartitionedCall|
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:??????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:?????????d: : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
f
J__inference_up_sampling2d_layer_call_and_return_conditional_losses_5331318

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"      2
Constc
Const_1Const*
_output_shapes
:*
dtype0*
valueB"      2	
Const_1X
mulMulConst:output:0Const_1:output:0*
T0*
_output_shapes
:2
mul?
resize/ResizeNearestNeighborResizeNearestNeighborinputsmul:z:0*
T0*/
_output_shapes
:?????????@*
half_pixel_centers(2
resize/ResizeNearestNeighbor?
IdentityIdentity-resize/ResizeNearestNeighbor:resized_images:0*
T0*/
_output_shapes
:?????????@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????@:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
?
b
F__inference_flatten_1_layer_call_and_return_conditional_losses_5329845

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"????    2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:??????????@2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:??????????@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:??????????:X T
0
_output_shapes
:??????????
 
_user_specified_nameinputs
?
G
+__inference_dropout_1_layer_call_fn_5331150

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_dropout_1_layer_call_and_return_conditional_losses_53296822
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:??????????:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
K
/__inference_up_sampling2d_layer_call_fn_5331328

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *S
fNRL
J__inference_up_sampling2d_layer_call_and_return_conditional_losses_53297452
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:?????????@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????@:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5329829

inputs&
readvariableop_resource:	?(
readvariableop_1_resource:	?7
(fusedbatchnormv3_readvariableop_resource:	?9
*fusedbatchnormv3_readvariableop_1_resource:	?
identity??FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1u
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp{
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes	
:?*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes	
:?*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*P
_output_shapes>
<:??????????:?:?:?:?:*
epsilon%o?:*
is_training( 2
FusedBatchNormV3x
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*0
_output_shapes
:??????????2

Identity?
NoOpNoOp ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:??????????: : : : 2B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:X T
0
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5331580

inputs&
readvariableop_resource:	?(
readvariableop_1_resource:	?7
(fusedbatchnormv3_readvariableop_resource:	?9
*fusedbatchnormv3_readvariableop_1_resource:	?
identity??FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1u
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp{
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes	
:?*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes	
:?*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*P
_output_shapes>
<:??????????:?:?:?:?:*
epsilon%o?:*
is_training( 2
FusedBatchNormV3x
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*0
_output_shapes
:??????????2

Identity?
NoOpNoOp ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:??????????: : : : 2B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:X T
0
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5331400

inputs%
readvariableop_resource:@'
readvariableop_1_resource:@6
(fusedbatchnormv3_readvariableop_resource:@8
*fusedbatchnormv3_readvariableop_1_resource:@
identity??AssignNewValue?AssignNewValue_1?FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1t
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
:@*
dtype02
ReadVariableOpz
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes
:@*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes
:@*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes
:@*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*K
_output_shapes9
7:?????????@:@:@:@:@:*
epsilon%o?:*
exponential_avg_factor%
?#<2
FusedBatchNormV3?
AssignNewValueAssignVariableOp(fusedbatchnormv3_readvariableop_resourceFusedBatchNormV3:batch_mean:0 ^FusedBatchNormV3/ReadVariableOp*
_output_shapes
 *
dtype02
AssignNewValue?
AssignNewValue_1AssignVariableOp*fusedbatchnormv3_readvariableop_1_resource!FusedBatchNormV3:batch_variance:0"^FusedBatchNormV3/ReadVariableOp_1*
_output_shapes
 *
dtype02
AssignNewValue_1w
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*/
_output_shapes
:?????????@2

Identity?
NoOpNoOp^AssignNewValue^AssignNewValue_1 ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*6
_input_shapes%
#:?????????@: : : : 2 
AssignNewValueAssignNewValue2$
AssignNewValue_1AssignNewValue_12B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
?
?
.__inference_sequential_1_layer_call_fn_5329968
dense_2_input
unknown:	d?
	unknown_0:	?
	unknown_1:	?
	unknown_2:	?
	unknown_3:	?
	unknown_4:	?#
	unknown_5:@@
	unknown_6:@
	unknown_7:@
	unknown_8:@
	unknown_9:@

unknown_10:@%

unknown_11:@?

unknown_12:	?

unknown_13:	?

unknown_14:	?

unknown_15:	?

unknown_16:	?

unknown_17:	?@2

unknown_18:2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCalldense_2_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18* 
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????2*6
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_sequential_1_layer_call_and_return_conditional_losses_53299252
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????22

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*N
_input_shapes=
;:?????????d: : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:V R
'
_output_shapes
:?????????d
'
_user_specified_namedense_2_input
?
G
+__inference_flatten_1_layer_call_fn_5331661

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_flatten_1_layer_call_and_return_conditional_losses_53298452
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:??????????@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:??????????:X T
0
_output_shapes
:??????????
 
_user_specified_nameinputs
??
?
"__inference__wrapped_model_5329155
dense_2_inputF
3sequential_1_dense_2_matmul_readvariableop_resource:	d?C
4sequential_1_dense_2_biasadd_readvariableop_resource:	?S
Dsequential_1_batch_normalization_3_batchnorm_readvariableop_resource:	?W
Hsequential_1_batch_normalization_3_batchnorm_mul_readvariableop_resource:	?U
Fsequential_1_batch_normalization_3_batchnorm_readvariableop_1_resource:	?U
Fsequential_1_batch_normalization_3_batchnorm_readvariableop_2_resource:	?N
4sequential_1_conv2d_2_conv2d_readvariableop_resource:@@C
5sequential_1_conv2d_2_biasadd_readvariableop_resource:@H
:sequential_1_batch_normalization_4_readvariableop_resource:@J
<sequential_1_batch_normalization_4_readvariableop_1_resource:@Y
Ksequential_1_batch_normalization_4_fusedbatchnormv3_readvariableop_resource:@[
Msequential_1_batch_normalization_4_fusedbatchnormv3_readvariableop_1_resource:@O
4sequential_1_conv2d_3_conv2d_readvariableop_resource:@?D
5sequential_1_conv2d_3_biasadd_readvariableop_resource:	?I
:sequential_1_batch_normalization_5_readvariableop_resource:	?K
<sequential_1_batch_normalization_5_readvariableop_1_resource:	?Z
Ksequential_1_batch_normalization_5_fusedbatchnormv3_readvariableop_resource:	?\
Msequential_1_batch_normalization_5_fusedbatchnormv3_readvariableop_1_resource:	?F
3sequential_1_dense_3_matmul_readvariableop_resource:	?@2B
4sequential_1_dense_3_biasadd_readvariableop_resource:2
identity??;sequential_1/batch_normalization_3/batchnorm/ReadVariableOp?=sequential_1/batch_normalization_3/batchnorm/ReadVariableOp_1?=sequential_1/batch_normalization_3/batchnorm/ReadVariableOp_2??sequential_1/batch_normalization_3/batchnorm/mul/ReadVariableOp?Bsequential_1/batch_normalization_4/FusedBatchNormV3/ReadVariableOp?Dsequential_1/batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1?1sequential_1/batch_normalization_4/ReadVariableOp?3sequential_1/batch_normalization_4/ReadVariableOp_1?Bsequential_1/batch_normalization_5/FusedBatchNormV3/ReadVariableOp?Dsequential_1/batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1?1sequential_1/batch_normalization_5/ReadVariableOp?3sequential_1/batch_normalization_5/ReadVariableOp_1?,sequential_1/conv2d_2/BiasAdd/ReadVariableOp?+sequential_1/conv2d_2/Conv2D/ReadVariableOp?,sequential_1/conv2d_3/BiasAdd/ReadVariableOp?+sequential_1/conv2d_3/Conv2D/ReadVariableOp?+sequential_1/dense_2/BiasAdd/ReadVariableOp?*sequential_1/dense_2/MatMul/ReadVariableOp?+sequential_1/dense_3/BiasAdd/ReadVariableOp?*sequential_1/dense_3/MatMul/ReadVariableOp?
*sequential_1/dense_2/MatMul/ReadVariableOpReadVariableOp3sequential_1_dense_2_matmul_readvariableop_resource*
_output_shapes
:	d?*
dtype02,
*sequential_1/dense_2/MatMul/ReadVariableOp?
sequential_1/dense_2/MatMulMatMuldense_2_input2sequential_1/dense_2/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
sequential_1/dense_2/MatMul?
+sequential_1/dense_2/BiasAdd/ReadVariableOpReadVariableOp4sequential_1_dense_2_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02-
+sequential_1/dense_2/BiasAdd/ReadVariableOp?
sequential_1/dense_2/BiasAddBiasAdd%sequential_1/dense_2/MatMul:product:03sequential_1/dense_2/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
sequential_1/dense_2/BiasAdd?
sequential_1/dense_2/EluElu%sequential_1/dense_2/BiasAdd:output:0*
T0*(
_output_shapes
:??????????2
sequential_1/dense_2/Elu?
sequential_1/dropout_1/IdentityIdentity&sequential_1/dense_2/Elu:activations:0*
T0*(
_output_shapes
:??????????2!
sequential_1/dropout_1/Identity?
;sequential_1/batch_normalization_3/batchnorm/ReadVariableOpReadVariableOpDsequential_1_batch_normalization_3_batchnorm_readvariableop_resource*
_output_shapes	
:?*
dtype02=
;sequential_1/batch_normalization_3/batchnorm/ReadVariableOp?
2sequential_1/batch_normalization_3/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:24
2sequential_1/batch_normalization_3/batchnorm/add/y?
0sequential_1/batch_normalization_3/batchnorm/addAddV2Csequential_1/batch_normalization_3/batchnorm/ReadVariableOp:value:0;sequential_1/batch_normalization_3/batchnorm/add/y:output:0*
T0*
_output_shapes	
:?22
0sequential_1/batch_normalization_3/batchnorm/add?
2sequential_1/batch_normalization_3/batchnorm/RsqrtRsqrt4sequential_1/batch_normalization_3/batchnorm/add:z:0*
T0*
_output_shapes	
:?24
2sequential_1/batch_normalization_3/batchnorm/Rsqrt?
?sequential_1/batch_normalization_3/batchnorm/mul/ReadVariableOpReadVariableOpHsequential_1_batch_normalization_3_batchnorm_mul_readvariableop_resource*
_output_shapes	
:?*
dtype02A
?sequential_1/batch_normalization_3/batchnorm/mul/ReadVariableOp?
0sequential_1/batch_normalization_3/batchnorm/mulMul6sequential_1/batch_normalization_3/batchnorm/Rsqrt:y:0Gsequential_1/batch_normalization_3/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:?22
0sequential_1/batch_normalization_3/batchnorm/mul?
2sequential_1/batch_normalization_3/batchnorm/mul_1Mul(sequential_1/dropout_1/Identity:output:04sequential_1/batch_normalization_3/batchnorm/mul:z:0*
T0*(
_output_shapes
:??????????24
2sequential_1/batch_normalization_3/batchnorm/mul_1?
=sequential_1/batch_normalization_3/batchnorm/ReadVariableOp_1ReadVariableOpFsequential_1_batch_normalization_3_batchnorm_readvariableop_1_resource*
_output_shapes	
:?*
dtype02?
=sequential_1/batch_normalization_3/batchnorm/ReadVariableOp_1?
2sequential_1/batch_normalization_3/batchnorm/mul_2MulEsequential_1/batch_normalization_3/batchnorm/ReadVariableOp_1:value:04sequential_1/batch_normalization_3/batchnorm/mul:z:0*
T0*
_output_shapes	
:?24
2sequential_1/batch_normalization_3/batchnorm/mul_2?
=sequential_1/batch_normalization_3/batchnorm/ReadVariableOp_2ReadVariableOpFsequential_1_batch_normalization_3_batchnorm_readvariableop_2_resource*
_output_shapes	
:?*
dtype02?
=sequential_1/batch_normalization_3/batchnorm/ReadVariableOp_2?
0sequential_1/batch_normalization_3/batchnorm/subSubEsequential_1/batch_normalization_3/batchnorm/ReadVariableOp_2:value:06sequential_1/batch_normalization_3/batchnorm/mul_2:z:0*
T0*
_output_shapes	
:?22
0sequential_1/batch_normalization_3/batchnorm/sub?
2sequential_1/batch_normalization_3/batchnorm/add_1AddV26sequential_1/batch_normalization_3/batchnorm/mul_1:z:04sequential_1/batch_normalization_3/batchnorm/sub:z:0*
T0*(
_output_shapes
:??????????24
2sequential_1/batch_normalization_3/batchnorm/add_1?
sequential_1/reshape_1/ShapeShape6sequential_1/batch_normalization_3/batchnorm/add_1:z:0*
T0*
_output_shapes
:2
sequential_1/reshape_1/Shape?
*sequential_1/reshape_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2,
*sequential_1/reshape_1/strided_slice/stack?
,sequential_1/reshape_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2.
,sequential_1/reshape_1/strided_slice/stack_1?
,sequential_1/reshape_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2.
,sequential_1/reshape_1/strided_slice/stack_2?
$sequential_1/reshape_1/strided_sliceStridedSlice%sequential_1/reshape_1/Shape:output:03sequential_1/reshape_1/strided_slice/stack:output:05sequential_1/reshape_1/strided_slice/stack_1:output:05sequential_1/reshape_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2&
$sequential_1/reshape_1/strided_slice?
&sequential_1/reshape_1/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2(
&sequential_1/reshape_1/Reshape/shape/1?
&sequential_1/reshape_1/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :2(
&sequential_1/reshape_1/Reshape/shape/2?
&sequential_1/reshape_1/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :@2(
&sequential_1/reshape_1/Reshape/shape/3?
$sequential_1/reshape_1/Reshape/shapePack-sequential_1/reshape_1/strided_slice:output:0/sequential_1/reshape_1/Reshape/shape/1:output:0/sequential_1/reshape_1/Reshape/shape/2:output:0/sequential_1/reshape_1/Reshape/shape/3:output:0*
N*
T0*
_output_shapes
:2&
$sequential_1/reshape_1/Reshape/shape?
sequential_1/reshape_1/ReshapeReshape6sequential_1/batch_normalization_3/batchnorm/add_1:z:0-sequential_1/reshape_1/Reshape/shape:output:0*
T0*/
_output_shapes
:?????????@2 
sequential_1/reshape_1/Reshape?
+sequential_1/conv2d_2/Conv2D/ReadVariableOpReadVariableOp4sequential_1_conv2d_2_conv2d_readvariableop_resource*&
_output_shapes
:@@*
dtype02-
+sequential_1/conv2d_2/Conv2D/ReadVariableOp?
sequential_1/conv2d_2/Conv2DConv2D'sequential_1/reshape_1/Reshape:output:03sequential_1/conv2d_2/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????@*
paddingSAME*
strides
2
sequential_1/conv2d_2/Conv2D?
,sequential_1/conv2d_2/BiasAdd/ReadVariableOpReadVariableOp5sequential_1_conv2d_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02.
,sequential_1/conv2d_2/BiasAdd/ReadVariableOp?
sequential_1/conv2d_2/BiasAddBiasAdd%sequential_1/conv2d_2/Conv2D:output:04sequential_1/conv2d_2/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????@2
sequential_1/conv2d_2/BiasAdd?
sequential_1/conv2d_2/EluElu&sequential_1/conv2d_2/BiasAdd:output:0*
T0*/
_output_shapes
:?????????@2
sequential_1/conv2d_2/Elu?
 sequential_1/up_sampling2d/ConstConst*
_output_shapes
:*
dtype0*
valueB"      2"
 sequential_1/up_sampling2d/Const?
"sequential_1/up_sampling2d/Const_1Const*
_output_shapes
:*
dtype0*
valueB"      2$
"sequential_1/up_sampling2d/Const_1?
sequential_1/up_sampling2d/mulMul)sequential_1/up_sampling2d/Const:output:0+sequential_1/up_sampling2d/Const_1:output:0*
T0*
_output_shapes
:2 
sequential_1/up_sampling2d/mul?
7sequential_1/up_sampling2d/resize/ResizeNearestNeighborResizeNearestNeighbor'sequential_1/conv2d_2/Elu:activations:0"sequential_1/up_sampling2d/mul:z:0*
T0*/
_output_shapes
:?????????@*
half_pixel_centers(29
7sequential_1/up_sampling2d/resize/ResizeNearestNeighbor?
1sequential_1/batch_normalization_4/ReadVariableOpReadVariableOp:sequential_1_batch_normalization_4_readvariableop_resource*
_output_shapes
:@*
dtype023
1sequential_1/batch_normalization_4/ReadVariableOp?
3sequential_1/batch_normalization_4/ReadVariableOp_1ReadVariableOp<sequential_1_batch_normalization_4_readvariableop_1_resource*
_output_shapes
:@*
dtype025
3sequential_1/batch_normalization_4/ReadVariableOp_1?
Bsequential_1/batch_normalization_4/FusedBatchNormV3/ReadVariableOpReadVariableOpKsequential_1_batch_normalization_4_fusedbatchnormv3_readvariableop_resource*
_output_shapes
:@*
dtype02D
Bsequential_1/batch_normalization_4/FusedBatchNormV3/ReadVariableOp?
Dsequential_1/batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1ReadVariableOpMsequential_1_batch_normalization_4_fusedbatchnormv3_readvariableop_1_resource*
_output_shapes
:@*
dtype02F
Dsequential_1/batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1?
3sequential_1/batch_normalization_4/FusedBatchNormV3FusedBatchNormV3Hsequential_1/up_sampling2d/resize/ResizeNearestNeighbor:resized_images:09sequential_1/batch_normalization_4/ReadVariableOp:value:0;sequential_1/batch_normalization_4/ReadVariableOp_1:value:0Jsequential_1/batch_normalization_4/FusedBatchNormV3/ReadVariableOp:value:0Lsequential_1/batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*K
_output_shapes9
7:?????????@:@:@:@:@:*
epsilon%o?:*
is_training( 25
3sequential_1/batch_normalization_4/FusedBatchNormV3?
+sequential_1/conv2d_3/Conv2D/ReadVariableOpReadVariableOp4sequential_1_conv2d_3_conv2d_readvariableop_resource*'
_output_shapes
:@?*
dtype02-
+sequential_1/conv2d_3/Conv2D/ReadVariableOp?
sequential_1/conv2d_3/Conv2DConv2D7sequential_1/batch_normalization_4/FusedBatchNormV3:y:03sequential_1/conv2d_3/Conv2D/ReadVariableOp:value:0*
T0*0
_output_shapes
:??????????*
paddingSAME*
strides
2
sequential_1/conv2d_3/Conv2D?
,sequential_1/conv2d_3/BiasAdd/ReadVariableOpReadVariableOp5sequential_1_conv2d_3_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02.
,sequential_1/conv2d_3/BiasAdd/ReadVariableOp?
sequential_1/conv2d_3/BiasAddBiasAdd%sequential_1/conv2d_3/Conv2D:output:04sequential_1/conv2d_3/BiasAdd/ReadVariableOp:value:0*
T0*0
_output_shapes
:??????????2
sequential_1/conv2d_3/BiasAdd?
sequential_1/conv2d_3/EluElu&sequential_1/conv2d_3/BiasAdd:output:0*
T0*0
_output_shapes
:??????????2
sequential_1/conv2d_3/Elu?
"sequential_1/up_sampling2d_1/ConstConst*
_output_shapes
:*
dtype0*
valueB"      2$
"sequential_1/up_sampling2d_1/Const?
$sequential_1/up_sampling2d_1/Const_1Const*
_output_shapes
:*
dtype0*
valueB"      2&
$sequential_1/up_sampling2d_1/Const_1?
 sequential_1/up_sampling2d_1/mulMul+sequential_1/up_sampling2d_1/Const:output:0-sequential_1/up_sampling2d_1/Const_1:output:0*
T0*
_output_shapes
:2"
 sequential_1/up_sampling2d_1/mul?
9sequential_1/up_sampling2d_1/resize/ResizeNearestNeighborResizeNearestNeighbor'sequential_1/conv2d_3/Elu:activations:0$sequential_1/up_sampling2d_1/mul:z:0*
T0*0
_output_shapes
:??????????*
half_pixel_centers(2;
9sequential_1/up_sampling2d_1/resize/ResizeNearestNeighbor?
1sequential_1/batch_normalization_5/ReadVariableOpReadVariableOp:sequential_1_batch_normalization_5_readvariableop_resource*
_output_shapes	
:?*
dtype023
1sequential_1/batch_normalization_5/ReadVariableOp?
3sequential_1/batch_normalization_5/ReadVariableOp_1ReadVariableOp<sequential_1_batch_normalization_5_readvariableop_1_resource*
_output_shapes	
:?*
dtype025
3sequential_1/batch_normalization_5/ReadVariableOp_1?
Bsequential_1/batch_normalization_5/FusedBatchNormV3/ReadVariableOpReadVariableOpKsequential_1_batch_normalization_5_fusedbatchnormv3_readvariableop_resource*
_output_shapes	
:?*
dtype02D
Bsequential_1/batch_normalization_5/FusedBatchNormV3/ReadVariableOp?
Dsequential_1/batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1ReadVariableOpMsequential_1_batch_normalization_5_fusedbatchnormv3_readvariableop_1_resource*
_output_shapes	
:?*
dtype02F
Dsequential_1/batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1?
3sequential_1/batch_normalization_5/FusedBatchNormV3FusedBatchNormV3Jsequential_1/up_sampling2d_1/resize/ResizeNearestNeighbor:resized_images:09sequential_1/batch_normalization_5/ReadVariableOp:value:0;sequential_1/batch_normalization_5/ReadVariableOp_1:value:0Jsequential_1/batch_normalization_5/FusedBatchNormV3/ReadVariableOp:value:0Lsequential_1/batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*P
_output_shapes>
<:??????????:?:?:?:?:*
epsilon%o?:*
is_training( 25
3sequential_1/batch_normalization_5/FusedBatchNormV3?
sequential_1/flatten_1/ConstConst*
_output_shapes
:*
dtype0*
valueB"????    2
sequential_1/flatten_1/Const?
sequential_1/flatten_1/ReshapeReshape7sequential_1/batch_normalization_5/FusedBatchNormV3:y:0%sequential_1/flatten_1/Const:output:0*
T0*(
_output_shapes
:??????????@2 
sequential_1/flatten_1/Reshape?
*sequential_1/dense_3/MatMul/ReadVariableOpReadVariableOp3sequential_1_dense_3_matmul_readvariableop_resource*
_output_shapes
:	?@2*
dtype02,
*sequential_1/dense_3/MatMul/ReadVariableOp?
sequential_1/dense_3/MatMulMatMul'sequential_1/flatten_1/Reshape:output:02sequential_1/dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????22
sequential_1/dense_3/MatMul?
+sequential_1/dense_3/BiasAdd/ReadVariableOpReadVariableOp4sequential_1_dense_3_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype02-
+sequential_1/dense_3/BiasAdd/ReadVariableOp?
sequential_1/dense_3/BiasAddBiasAdd%sequential_1/dense_3/MatMul:product:03sequential_1/dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????22
sequential_1/dense_3/BiasAdd?
sequential_1/dense_3/SigmoidSigmoid%sequential_1/dense_3/BiasAdd:output:0*
T0*'
_output_shapes
:?????????22
sequential_1/dense_3/Sigmoid{
IdentityIdentity sequential_1/dense_3/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:?????????22

Identity?	
NoOpNoOp<^sequential_1/batch_normalization_3/batchnorm/ReadVariableOp>^sequential_1/batch_normalization_3/batchnorm/ReadVariableOp_1>^sequential_1/batch_normalization_3/batchnorm/ReadVariableOp_2@^sequential_1/batch_normalization_3/batchnorm/mul/ReadVariableOpC^sequential_1/batch_normalization_4/FusedBatchNormV3/ReadVariableOpE^sequential_1/batch_normalization_4/FusedBatchNormV3/ReadVariableOp_12^sequential_1/batch_normalization_4/ReadVariableOp4^sequential_1/batch_normalization_4/ReadVariableOp_1C^sequential_1/batch_normalization_5/FusedBatchNormV3/ReadVariableOpE^sequential_1/batch_normalization_5/FusedBatchNormV3/ReadVariableOp_12^sequential_1/batch_normalization_5/ReadVariableOp4^sequential_1/batch_normalization_5/ReadVariableOp_1-^sequential_1/conv2d_2/BiasAdd/ReadVariableOp,^sequential_1/conv2d_2/Conv2D/ReadVariableOp-^sequential_1/conv2d_3/BiasAdd/ReadVariableOp,^sequential_1/conv2d_3/Conv2D/ReadVariableOp,^sequential_1/dense_2/BiasAdd/ReadVariableOp+^sequential_1/dense_2/MatMul/ReadVariableOp,^sequential_1/dense_3/BiasAdd/ReadVariableOp+^sequential_1/dense_3/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*N
_input_shapes=
;:?????????d: : : : : : : : : : : : : : : : : : : : 2z
;sequential_1/batch_normalization_3/batchnorm/ReadVariableOp;sequential_1/batch_normalization_3/batchnorm/ReadVariableOp2~
=sequential_1/batch_normalization_3/batchnorm/ReadVariableOp_1=sequential_1/batch_normalization_3/batchnorm/ReadVariableOp_12~
=sequential_1/batch_normalization_3/batchnorm/ReadVariableOp_2=sequential_1/batch_normalization_3/batchnorm/ReadVariableOp_22?
?sequential_1/batch_normalization_3/batchnorm/mul/ReadVariableOp?sequential_1/batch_normalization_3/batchnorm/mul/ReadVariableOp2?
Bsequential_1/batch_normalization_4/FusedBatchNormV3/ReadVariableOpBsequential_1/batch_normalization_4/FusedBatchNormV3/ReadVariableOp2?
Dsequential_1/batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1Dsequential_1/batch_normalization_4/FusedBatchNormV3/ReadVariableOp_12f
1sequential_1/batch_normalization_4/ReadVariableOp1sequential_1/batch_normalization_4/ReadVariableOp2j
3sequential_1/batch_normalization_4/ReadVariableOp_13sequential_1/batch_normalization_4/ReadVariableOp_12?
Bsequential_1/batch_normalization_5/FusedBatchNormV3/ReadVariableOpBsequential_1/batch_normalization_5/FusedBatchNormV3/ReadVariableOp2?
Dsequential_1/batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1Dsequential_1/batch_normalization_5/FusedBatchNormV3/ReadVariableOp_12f
1sequential_1/batch_normalization_5/ReadVariableOp1sequential_1/batch_normalization_5/ReadVariableOp2j
3sequential_1/batch_normalization_5/ReadVariableOp_13sequential_1/batch_normalization_5/ReadVariableOp_12\
,sequential_1/conv2d_2/BiasAdd/ReadVariableOp,sequential_1/conv2d_2/BiasAdd/ReadVariableOp2Z
+sequential_1/conv2d_2/Conv2D/ReadVariableOp+sequential_1/conv2d_2/Conv2D/ReadVariableOp2\
,sequential_1/conv2d_3/BiasAdd/ReadVariableOp,sequential_1/conv2d_3/BiasAdd/ReadVariableOp2Z
+sequential_1/conv2d_3/Conv2D/ReadVariableOp+sequential_1/conv2d_3/Conv2D/ReadVariableOp2Z
+sequential_1/dense_2/BiasAdd/ReadVariableOp+sequential_1/dense_2/BiasAdd/ReadVariableOp2X
*sequential_1/dense_2/MatMul/ReadVariableOp*sequential_1/dense_2/MatMul/ReadVariableOp2Z
+sequential_1/dense_3/BiasAdd/ReadVariableOp+sequential_1/dense_3/BiasAdd/ReadVariableOp2X
*sequential_1/dense_3/MatMul/ReadVariableOp*sequential_1/dense_3/MatMul/ReadVariableOp:V R
'
_output_shapes
:?????????d
'
_user_specified_namedense_2_input
?
?
E__inference_conv2d_3_layer_call_and_return_conditional_losses_5329797

inputs9
conv2d_readvariableop_resource:@?.
biasadd_readvariableop_resource:	?
identity??BiasAdd/ReadVariableOp?Conv2D/ReadVariableOp?/conv2d_3/bias/Regularizer/Square/ReadVariableOp?1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*'
_output_shapes
:@?*
dtype02
Conv2D/ReadVariableOp?
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*0
_output_shapes
:??????????*
paddingSAME*
strides
2
Conv2D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*0
_output_shapes
:??????????2	
BiasAdd^
EluEluBiasAdd:output:0*
T0*0
_output_shapes
:??????????2
Elu?
1conv2d_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*'
_output_shapes
:@?*
dtype023
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_3/kernel/Regularizer/SquareSquare9conv2d_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*'
_output_shapes
:@?2$
"conv2d_3/kernel/Regularizer/Square?
!conv2d_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_3/kernel/Regularizer/Const?
conv2d_3/kernel/Regularizer/SumSum&conv2d_3/kernel/Regularizer/Square:y:0*conv2d_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/Sum?
!conv2d_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_3/kernel/Regularizer/mul/x?
conv2d_3/kernel/Regularizer/mulMul*conv2d_3/kernel/Regularizer/mul/x:output:0(conv2d_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/mul?
/conv2d_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:?*
dtype021
/conv2d_3/bias/Regularizer/Square/ReadVariableOp?
 conv2d_3/bias/Regularizer/SquareSquare7conv2d_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2"
 conv2d_3/bias/Regularizer/Square?
conv2d_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_3/bias/Regularizer/Const?
conv2d_3/bias/Regularizer/SumSum$conv2d_3/bias/Regularizer/Square:y:0(conv2d_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/Sum?
conv2d_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_3/bias/Regularizer/mul/x?
conv2d_3/bias/Regularizer/mulMul(conv2d_3/bias/Regularizer/mul/x:output:0&conv2d_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/mulu
IdentityIdentityElu:activations:0^NoOp*
T0*0
_output_shapes
:??????????2

Identity?
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp0^conv2d_3/bias/Regularizer/Square/ReadVariableOp2^conv2d_3/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp2b
/conv2d_3/bias/Regularizer/Square/ReadVariableOp/conv2d_3/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp1conv2d_3/kernel/Regularizer/Square/ReadVariableOp:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
?
e
F__inference_dropout_1_layer_call_and_return_conditional_losses_5330128

inputs
identity?c
dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout/Constt
dropout/MulMulinputsdropout/Const:output:0*
T0*(
_output_shapes
:??????????2
dropout/MulT
dropout/ShapeShapeinputs*
T0*
_output_shapes
:2
dropout/Shape?
$dropout/random_uniform/RandomUniformRandomUniformdropout/Shape:output:0*
T0*(
_output_shapes
:??????????*
dtype02&
$dropout/random_uniform/RandomUniformu
dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2
dropout/GreaterEqual/y?
dropout/GreaterEqualGreaterEqual-dropout/random_uniform/RandomUniform:output:0dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:??????????2
dropout/GreaterEqual?
dropout/CastCastdropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:??????????2
dropout/Cast{
dropout/Mul_1Muldropout/Mul:z:0dropout/Cast:y:0*
T0*(
_output_shapes
:??????????2
dropout/Mul_1f
IdentityIdentitydropout/Mul_1:z:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:??????????:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5329419

inputs%
readvariableop_resource:@'
readvariableop_1_resource:@6
(fusedbatchnormv3_readvariableop_resource:@8
*fusedbatchnormv3_readvariableop_1_resource:@
identity??AssignNewValue?AssignNewValue_1?FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1t
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
:@*
dtype02
ReadVariableOpz
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes
:@*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes
:@*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes
:@*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*]
_output_shapesK
I:+???????????????????????????@:@:@:@:@:*
epsilon%o?:*
exponential_avg_factor%
?#<2
FusedBatchNormV3?
AssignNewValueAssignVariableOp(fusedbatchnormv3_readvariableop_resourceFusedBatchNormV3:batch_mean:0 ^FusedBatchNormV3/ReadVariableOp*
_output_shapes
 *
dtype02
AssignNewValue?
AssignNewValue_1AssignVariableOp*fusedbatchnormv3_readvariableop_1_resource!FusedBatchNormV3:batch_variance:0"^FusedBatchNormV3/ReadVariableOp_1*
_output_shapes
 *
dtype02
AssignNewValue_1?
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*A
_output_shapes/
-:+???????????????????????????@2

Identity?
NoOpNoOp^AssignNewValue^AssignNewValue_1 ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*H
_input_shapes7
5:+???????????????????????????@: : : : 2 
AssignNewValueAssignNewValue2$
AssignNewValue_1AssignNewValue_12B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:i e
A
_output_shapes/
-:+???????????????????????????@
 
_user_specified_nameinputs
?
?
7__inference_batch_normalization_4_layer_call_fn_5331439

inputs
unknown:@
	unknown_0:@
	unknown_1:@
	unknown_2:@
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_53297642
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:?????????@2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*6
_input_shapes%
#:?????????@: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
?
b
F__inference_reshape_1_layer_call_and_return_conditional_losses_5331249

inputs
identityD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2?
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_sliced
Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2
Reshape/shape/1d
Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :2
Reshape/shape/2d
Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :@2
Reshape/shape/3?
Reshape/shapePackstrided_slice:output:0Reshape/shape/1:output:0Reshape/shape/2:output:0Reshape/shape/3:output:0*
N*
T0*
_output_shapes
:2
Reshape/shapew
ReshapeReshapeinputsReshape/shape:output:0*
T0*/
_output_shapes
:?????????@2	
Reshapel
IdentityIdentityReshape:output:0*
T0*/
_output_shapes
:?????????@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:??????????:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
b
F__inference_flatten_1_layer_call_and_return_conditional_losses_5331656

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"????    2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:??????????@2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:??????????@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:??????????:X T
0
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5329581

inputs&
readvariableop_resource:	?(
readvariableop_1_resource:	?7
(fusedbatchnormv3_readvariableop_resource:	?9
*fusedbatchnormv3_readvariableop_1_resource:	?
identity??AssignNewValue?AssignNewValue_1?FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1u
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp{
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes	
:?*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes	
:?*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*b
_output_shapesP
N:,????????????????????????????:?:?:?:?:*
epsilon%o?:*
exponential_avg_factor%
?#<2
FusedBatchNormV3?
AssignNewValueAssignVariableOp(fusedbatchnormv3_readvariableop_resourceFusedBatchNormV3:batch_mean:0 ^FusedBatchNormV3/ReadVariableOp*
_output_shapes
 *
dtype02
AssignNewValue?
AssignNewValue_1AssignVariableOp*fusedbatchnormv3_readvariableop_1_resource!FusedBatchNormV3:batch_variance:0"^FusedBatchNormV3/ReadVariableOp_1*
_output_shapes
 *
dtype02
AssignNewValue_1?
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*B
_output_shapes0
.:,????????????????????????????2

Identity?
NoOpNoOp^AssignNewValue^AssignNewValue_1 ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:,????????????????????????????: : : : 2 
AssignNewValueAssignNewValue2$
AssignNewValue_1AssignNewValue_12B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:j f
B
_output_shapes0
.:,????????????????????????????
 
_user_specified_nameinputs
?
?
D__inference_dense_3_layer_call_and_return_conditional_losses_5329870

inputs1
matmul_readvariableop_resource:	?@2-
biasadd_readvariableop_resource:2
identity??BiasAdd/ReadVariableOp?MatMul/ReadVariableOp?.dense_3/bias/Regularizer/Square/ReadVariableOp?0dense_3/kernel/Regularizer/Square/ReadVariableOp?
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	?@2*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????22
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????22	
BiasAdda
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:?????????22	
Sigmoid?
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	?@2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp?
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?@22#
!dense_3/kernel/Regularizer/Square?
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const?
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum?
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_3/kernel/Regularizer/mul/x?
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul?
.dense_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype020
.dense_3/bias/Regularizer/Square/ReadVariableOp?
dense_3/bias/Regularizer/SquareSquare6dense_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:22!
dense_3/bias/Regularizer/Square?
dense_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_3/bias/Regularizer/Const?
dense_3/bias/Regularizer/SumSum#dense_3/bias/Regularizer/Square:y:0'dense_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/Sum?
dense_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_3/bias/Regularizer/mul/x?
dense_3/bias/Regularizer/mulMul'dense_3/bias/Regularizer/mul/x:output:0%dense_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/mulf
IdentityIdentitySigmoid:y:0^NoOp*
T0*'
_output_shapes
:?????????22

Identity?
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp/^dense_3/bias/Regularizer/Square/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:??????????@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2`
.dense_3/bias/Regularizer/Square/ReadVariableOp.dense_3/bias/Regularizer/Square/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:??????????@
 
_user_specified_nameinputs
?
h
L__inference_up_sampling2d_1_layer_call_and_return_conditional_losses_5329495

inputs
identityD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2?
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2
strided_slice_
ConstConst*
_output_shapes
:*
dtype0*
valueB"      2
Const^
mulMulstrided_slice:output:0Const:output:0*
T0*
_output_shapes
:2
mul?
resize/ResizeNearestNeighborResizeNearestNeighborinputsmul:z:0*
T0*J
_output_shapes8
6:4????????????????????????????????????*
half_pixel_centers(2
resize/ResizeNearestNeighbor?
IdentityIdentity-resize/ResizeNearestNeighbor:resized_images:0*
T0*J
_output_shapes8
6:4????????????????????????????????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4????????????????????????????????????:r n
J
_output_shapes8
6:4????????????????????????????????????
 
_user_specified_nameinputs
?
?
D__inference_dense_3_layer_call_and_return_conditional_losses_5331696

inputs1
matmul_readvariableop_resource:	?@2-
biasadd_readvariableop_resource:2
identity??BiasAdd/ReadVariableOp?MatMul/ReadVariableOp?.dense_3/bias/Regularizer/Square/ReadVariableOp?0dense_3/kernel/Regularizer/Square/ReadVariableOp?
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	?@2*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????22
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????22	
BiasAdda
SigmoidSigmoidBiasAdd:output:0*
T0*'
_output_shapes
:?????????22	
Sigmoid?
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	?@2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp?
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?@22#
!dense_3/kernel/Regularizer/Square?
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const?
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum?
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_3/kernel/Regularizer/mul/x?
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul?
.dense_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:2*
dtype020
.dense_3/bias/Regularizer/Square/ReadVariableOp?
dense_3/bias/Regularizer/SquareSquare6dense_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:22!
dense_3/bias/Regularizer/Square?
dense_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_3/bias/Regularizer/Const?
dense_3/bias/Regularizer/SumSum#dense_3/bias/Regularizer/Square:y:0'dense_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/Sum?
dense_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_3/bias/Regularizer/mul/x?
dense_3/bias/Regularizer/mulMul'dense_3/bias/Regularizer/mul/x:output:0%dense_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/mulf
IdentityIdentitySigmoid:y:0^NoOp*
T0*'
_output_shapes
:?????????22

Identity?
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp/^dense_3/bias/Regularizer/Square/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:??????????@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2`
.dense_3/bias/Regularizer/Square/ReadVariableOp.dense_3/bias/Regularizer/Square/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp:P L
(
_output_shapes
:??????????@
 
_user_specified_nameinputs
?	
?
7__inference_batch_normalization_4_layer_call_fn_5331413

inputs
unknown:@
	unknown_0:@
	unknown_1:@
	unknown_2:@
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *A
_output_shapes/
-:+???????????????????????????@*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_53293752
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*A
_output_shapes/
-:+???????????????????????????@2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*H
_input_shapes7
5:+???????????????????????????@: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:i e
A
_output_shapes/
-:+???????????????????????????@
 
_user_specified_nameinputs
?
?
E__inference_conv2d_3_layer_call_and_return_conditional_losses_5331487

inputs9
conv2d_readvariableop_resource:@?.
biasadd_readvariableop_resource:	?
identity??BiasAdd/ReadVariableOp?Conv2D/ReadVariableOp?/conv2d_3/bias/Regularizer/Square/ReadVariableOp?1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*'
_output_shapes
:@?*
dtype02
Conv2D/ReadVariableOp?
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*0
_output_shapes
:??????????*
paddingSAME*
strides
2
Conv2D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*0
_output_shapes
:??????????2	
BiasAdd^
EluEluBiasAdd:output:0*
T0*0
_output_shapes
:??????????2
Elu?
1conv2d_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*'
_output_shapes
:@?*
dtype023
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_3/kernel/Regularizer/SquareSquare9conv2d_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*'
_output_shapes
:@?2$
"conv2d_3/kernel/Regularizer/Square?
!conv2d_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_3/kernel/Regularizer/Const?
conv2d_3/kernel/Regularizer/SumSum&conv2d_3/kernel/Regularizer/Square:y:0*conv2d_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/Sum?
!conv2d_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_3/kernel/Regularizer/mul/x?
conv2d_3/kernel/Regularizer/mulMul*conv2d_3/kernel/Regularizer/mul/x:output:0(conv2d_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/mul?
/conv2d_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:?*
dtype021
/conv2d_3/bias/Regularizer/Square/ReadVariableOp?
 conv2d_3/bias/Regularizer/SquareSquare7conv2d_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2"
 conv2d_3/bias/Regularizer/Square?
conv2d_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_3/bias/Regularizer/Const?
conv2d_3/bias/Regularizer/SumSum$conv2d_3/bias/Regularizer/Square:y:0(conv2d_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/Sum?
conv2d_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_3/bias/Regularizer/mul/x?
conv2d_3/bias/Regularizer/mulMul(conv2d_3/bias/Regularizer/mul/x:output:0&conv2d_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/mulu
IdentityIdentityElu:activations:0^NoOp*
T0*0
_output_shapes
:??????????2

Identity?
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp0^conv2d_3/bias/Regularizer/Square/ReadVariableOp2^conv2d_3/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp2b
/conv2d_3/bias/Regularizer/Square/ReadVariableOp/conv2d_3/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp1conv2d_3/kernel/Regularizer/Square/ReadVariableOp:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
?
M
1__inference_up_sampling2d_1_layer_call_fn_5331521

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *J
_output_shapes8
6:4????????????????????????????????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_up_sampling2d_1_layer_call_and_return_conditional_losses_53294952
PartitionedCall?
IdentityIdentityPartitionedCall:output:0*
T0*J
_output_shapes8
6:4????????????????????????????????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4????????????????????????????????????:r n
J
_output_shapes8
6:4????????????????????????????????????
 
_user_specified_nameinputs
?4
?	
 __inference__traced_save_5331876
file_prefix-
)savev2_dense_2_kernel_read_readvariableop+
'savev2_dense_2_bias_read_readvariableop:
6savev2_batch_normalization_3_gamma_read_readvariableop9
5savev2_batch_normalization_3_beta_read_readvariableop@
<savev2_batch_normalization_3_moving_mean_read_readvariableopD
@savev2_batch_normalization_3_moving_variance_read_readvariableop.
*savev2_conv2d_2_kernel_read_readvariableop,
(savev2_conv2d_2_bias_read_readvariableop:
6savev2_batch_normalization_4_gamma_read_readvariableop9
5savev2_batch_normalization_4_beta_read_readvariableop@
<savev2_batch_normalization_4_moving_mean_read_readvariableopD
@savev2_batch_normalization_4_moving_variance_read_readvariableop.
*savev2_conv2d_3_kernel_read_readvariableop,
(savev2_conv2d_3_bias_read_readvariableop:
6savev2_batch_normalization_5_gamma_read_readvariableop9
5savev2_batch_normalization_5_beta_read_readvariableop@
<savev2_batch_normalization_5_moving_mean_read_readvariableopD
@savev2_batch_normalization_5_moving_variance_read_readvariableop-
)savev2_dense_3_kernel_read_readvariableop+
'savev2_dense_3_bias_read_readvariableop
savev2_const

identity_1??MergeV2Checkpoints?
StaticRegexFullMatchStaticRegexFullMatchfile_prefix"/device:CPU:**
_output_shapes
: *
pattern
^s3://.*2
StaticRegexFullMatchc
ConstConst"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B.part2
Constl
Const_1Const"/device:CPU:**
_output_shapes
: *
dtype0*
valueB B
_temp/part2	
Const_1?
SelectSelectStaticRegexFullMatch:output:0Const:output:0Const_1:output:0"/device:CPU:**
T0*
_output_shapes
: 2
Selectt

StringJoin
StringJoinfile_prefixSelect:output:0"/device:CPU:**
N*
_output_shapes
: 2

StringJoinZ

num_shardsConst*
_output_shapes
: *
dtype0*
value	B :2

num_shards
ShardedFilename/shardConst"/device:CPU:0*
_output_shapes
: *
dtype0*
value	B : 2
ShardedFilename/shard?
ShardedFilenameShardedFilenameStringJoin:output:0ShardedFilename/shard:output:0num_shards:output:0"/device:CPU:0*
_output_shapes
: 2
ShardedFilename?

SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*?	
value?	B?	B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-1/gamma/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/beta/.ATTRIBUTES/VARIABLE_VALUEB;layer_with_weights-1/moving_mean/.ATTRIBUTES/VARIABLE_VALUEB?layer_with_weights-1/moving_variance/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-3/gamma/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/beta/.ATTRIBUTES/VARIABLE_VALUEB;layer_with_weights-3/moving_mean/.ATTRIBUTES/VARIABLE_VALUEB?layer_with_weights-3/moving_variance/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-5/gamma/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/beta/.ATTRIBUTES/VARIABLE_VALUEB;layer_with_weights-5/moving_mean/.ATTRIBUTES/VARIABLE_VALUEB?layer_with_weights-5/moving_variance/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
SaveV2/tensor_names?
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*=
value4B2B B B B B B B B B B B B B B B B B B B B B 2
SaveV2/shape_and_slices?	
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0)savev2_dense_2_kernel_read_readvariableop'savev2_dense_2_bias_read_readvariableop6savev2_batch_normalization_3_gamma_read_readvariableop5savev2_batch_normalization_3_beta_read_readvariableop<savev2_batch_normalization_3_moving_mean_read_readvariableop@savev2_batch_normalization_3_moving_variance_read_readvariableop*savev2_conv2d_2_kernel_read_readvariableop(savev2_conv2d_2_bias_read_readvariableop6savev2_batch_normalization_4_gamma_read_readvariableop5savev2_batch_normalization_4_beta_read_readvariableop<savev2_batch_normalization_4_moving_mean_read_readvariableop@savev2_batch_normalization_4_moving_variance_read_readvariableop*savev2_conv2d_3_kernel_read_readvariableop(savev2_conv2d_3_bias_read_readvariableop6savev2_batch_normalization_5_gamma_read_readvariableop5savev2_batch_normalization_5_beta_read_readvariableop<savev2_batch_normalization_5_moving_mean_read_readvariableop@savev2_batch_normalization_5_moving_variance_read_readvariableop)savev2_dense_3_kernel_read_readvariableop'savev2_dense_3_bias_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *#
dtypes
22
SaveV2?
&MergeV2Checkpoints/checkpoint_prefixesPackShardedFilename:filename:0^SaveV2"/device:CPU:0*
N*
T0*
_output_shapes
:2(
&MergeV2Checkpoints/checkpoint_prefixes?
MergeV2CheckpointsMergeV2Checkpoints/MergeV2Checkpoints/checkpoint_prefixes:output:0file_prefix"/device:CPU:0*
_output_shapes
 2
MergeV2Checkpointsr
IdentityIdentityfile_prefix^MergeV2Checkpoints"/device:CPU:0*
T0*
_output_shapes
: 2

Identity_

Identity_1IdentityIdentity:output:0^NoOp*
T0*
_output_shapes
: 2

Identity_1c
NoOpNoOp^MergeV2Checkpoints*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"!

identity_1Identity_1:output:0*?
_input_shapes?
?: :	d?:?:?:?:?:?:@@:@:@:@:@:@:@?:?:?:?:?:?:	?@2:2: 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:%!

_output_shapes
:	d?:!

_output_shapes	
:?:!

_output_shapes	
:?:!

_output_shapes	
:?:!

_output_shapes	
:?:!

_output_shapes	
:?:,(
&
_output_shapes
:@@: 

_output_shapes
:@: 	

_output_shapes
:@: 


_output_shapes
:@: 

_output_shapes
:@: 

_output_shapes
:@:-)
'
_output_shapes
:@?:!

_output_shapes	
:?:!

_output_shapes	
:?:!

_output_shapes	
:?:!

_output_shapes	
:?:!

_output_shapes	
:?:%!

_output_shapes
:	?@2: 

_output_shapes
:2:

_output_shapes
: 
??
?
I__inference_sequential_1_layer_call_and_return_conditional_losses_5330484
dense_2_input"
dense_2_5330383:	d?
dense_2_5330385:	?,
batch_normalization_3_5330389:	?,
batch_normalization_3_5330391:	?,
batch_normalization_3_5330393:	?,
batch_normalization_3_5330395:	?*
conv2d_2_5330399:@@
conv2d_2_5330401:@+
batch_normalization_4_5330405:@+
batch_normalization_4_5330407:@+
batch_normalization_4_5330409:@+
batch_normalization_4_5330411:@+
conv2d_3_5330414:@?
conv2d_3_5330416:	?,
batch_normalization_5_5330420:	?,
batch_normalization_5_5330422:	?,
batch_normalization_5_5330424:	?,
batch_normalization_5_5330426:	?"
dense_3_5330430:	?@2
dense_3_5330432:2
identity??-batch_normalization_3/StatefulPartitionedCall?-batch_normalization_4/StatefulPartitionedCall?-batch_normalization_5/StatefulPartitionedCall? conv2d_2/StatefulPartitionedCall?/conv2d_2/bias/Regularizer/Square/ReadVariableOp?1conv2d_2/kernel/Regularizer/Square/ReadVariableOp? conv2d_3/StatefulPartitionedCall?/conv2d_3/bias/Regularizer/Square/ReadVariableOp?1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?dense_2/StatefulPartitionedCall?.dense_2/bias/Regularizer/Square/ReadVariableOp?0dense_2/kernel/Regularizer/Square/ReadVariableOp?dense_3/StatefulPartitionedCall?.dense_3/bias/Regularizer/Square/ReadVariableOp?0dense_3/kernel/Regularizer/Square/ReadVariableOp?
dense_2/StatefulPartitionedCallStatefulPartitionedCalldense_2_inputdense_2_5330383dense_2_5330385*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_2_layer_call_and_return_conditional_losses_53296712!
dense_2/StatefulPartitionedCall?
dropout_1/PartitionedCallPartitionedCall(dense_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_dropout_1_layer_call_and_return_conditional_losses_53296822
dropout_1/PartitionedCall?
-batch_normalization_3/StatefulPartitionedCallStatefulPartitionedCall"dropout_1/PartitionedCall:output:0batch_normalization_3_5330389batch_normalization_3_5330391batch_normalization_3_5330393batch_normalization_3_5330395*
Tin	
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_53291792/
-batch_normalization_3/StatefulPartitionedCall?
reshape_1/PartitionedCallPartitionedCall6batch_normalization_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_reshape_1_layer_call_and_return_conditional_losses_53297072
reshape_1/PartitionedCall?
 conv2d_2/StatefulPartitionedCallStatefulPartitionedCall"reshape_1/PartitionedCall:output:0conv2d_2_5330399conv2d_2_5330401*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv2d_2_layer_call_and_return_conditional_losses_53297322"
 conv2d_2/StatefulPartitionedCall?
up_sampling2d/PartitionedCallPartitionedCall)conv2d_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *S
fNRL
J__inference_up_sampling2d_layer_call_and_return_conditional_losses_53297452
up_sampling2d/PartitionedCall?
-batch_normalization_4/StatefulPartitionedCallStatefulPartitionedCall&up_sampling2d/PartitionedCall:output:0batch_normalization_4_5330405batch_normalization_4_5330407batch_normalization_4_5330409batch_normalization_4_5330411*
Tin	
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_53297642/
-batch_normalization_4/StatefulPartitionedCall?
 conv2d_3/StatefulPartitionedCallStatefulPartitionedCall6batch_normalization_4/StatefulPartitionedCall:output:0conv2d_3_5330414conv2d_3_5330416*
Tin
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv2d_3_layer_call_and_return_conditional_losses_53297972"
 conv2d_3/StatefulPartitionedCall?
up_sampling2d_1/PartitionedCallPartitionedCall)conv2d_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_up_sampling2d_1_layer_call_and_return_conditional_losses_53298102!
up_sampling2d_1/PartitionedCall?
-batch_normalization_5/StatefulPartitionedCallStatefulPartitionedCall(up_sampling2d_1/PartitionedCall:output:0batch_normalization_5_5330420batch_normalization_5_5330422batch_normalization_5_5330424batch_normalization_5_5330426*
Tin	
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_53298292/
-batch_normalization_5/StatefulPartitionedCall?
flatten_1/PartitionedCallPartitionedCall6batch_normalization_5/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_flatten_1_layer_call_and_return_conditional_losses_53298452
flatten_1/PartitionedCall?
dense_3/StatefulPartitionedCallStatefulPartitionedCall"flatten_1/PartitionedCall:output:0dense_3_5330430dense_3_5330432*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????2*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_3_layer_call_and_return_conditional_losses_53298702!
dense_3/StatefulPartitionedCall?
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_5330383*
_output_shapes
:	d?*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp?
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	d?2#
!dense_2/kernel/Regularizer/Square?
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const?
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/Sum?
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_2/kernel/Regularizer/mul/x?
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mul?
.dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_5330385*
_output_shapes	
:?*
dtype020
.dense_2/bias/Regularizer/Square/ReadVariableOp?
dense_2/bias/Regularizer/SquareSquare6dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2!
dense_2/bias/Regularizer/Square?
dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_2/bias/Regularizer/Const?
dense_2/bias/Regularizer/SumSum#dense_2/bias/Regularizer/Square:y:0'dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/Sum?
dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_2/bias/Regularizer/mul/x?
dense_2/bias/Regularizer/mulMul'dense_2/bias/Regularizer/mul/x:output:0%dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/mul?
1conv2d_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_2_5330399*&
_output_shapes
:@@*
dtype023
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_2/kernel/Regularizer/SquareSquare9conv2d_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*&
_output_shapes
:@@2$
"conv2d_2/kernel/Regularizer/Square?
!conv2d_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_2/kernel/Regularizer/Const?
conv2d_2/kernel/Regularizer/SumSum&conv2d_2/kernel/Regularizer/Square:y:0*conv2d_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/Sum?
!conv2d_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_2/kernel/Regularizer/mul/x?
conv2d_2/kernel/Regularizer/mulMul*conv2d_2/kernel/Regularizer/mul/x:output:0(conv2d_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/mul?
/conv2d_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_2_5330401*
_output_shapes
:@*
dtype021
/conv2d_2/bias/Regularizer/Square/ReadVariableOp?
 conv2d_2/bias/Regularizer/SquareSquare7conv2d_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2"
 conv2d_2/bias/Regularizer/Square?
conv2d_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_2/bias/Regularizer/Const?
conv2d_2/bias/Regularizer/SumSum$conv2d_2/bias/Regularizer/Square:y:0(conv2d_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/Sum?
conv2d_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_2/bias/Regularizer/mul/x?
conv2d_2/bias/Regularizer/mulMul(conv2d_2/bias/Regularizer/mul/x:output:0&conv2d_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/mul?
1conv2d_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_3_5330414*'
_output_shapes
:@?*
dtype023
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_3/kernel/Regularizer/SquareSquare9conv2d_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*'
_output_shapes
:@?2$
"conv2d_3/kernel/Regularizer/Square?
!conv2d_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_3/kernel/Regularizer/Const?
conv2d_3/kernel/Regularizer/SumSum&conv2d_3/kernel/Regularizer/Square:y:0*conv2d_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/Sum?
!conv2d_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_3/kernel/Regularizer/mul/x?
conv2d_3/kernel/Regularizer/mulMul*conv2d_3/kernel/Regularizer/mul/x:output:0(conv2d_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/mul?
/conv2d_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_3_5330416*
_output_shapes	
:?*
dtype021
/conv2d_3/bias/Regularizer/Square/ReadVariableOp?
 conv2d_3/bias/Regularizer/SquareSquare7conv2d_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2"
 conv2d_3/bias/Regularizer/Square?
conv2d_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_3/bias/Regularizer/Const?
conv2d_3/bias/Regularizer/SumSum$conv2d_3/bias/Regularizer/Square:y:0(conv2d_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/Sum?
conv2d_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_3/bias/Regularizer/mul/x?
conv2d_3/bias/Regularizer/mulMul(conv2d_3/bias/Regularizer/mul/x:output:0&conv2d_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/mul?
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_5330430*
_output_shapes
:	?@2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp?
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?@22#
!dense_3/kernel/Regularizer/Square?
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const?
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum?
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_3/kernel/Regularizer/mul/x?
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul?
.dense_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_5330432*
_output_shapes
:2*
dtype020
.dense_3/bias/Regularizer/Square/ReadVariableOp?
dense_3/bias/Regularizer/SquareSquare6dense_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:22!
dense_3/bias/Regularizer/Square?
dense_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_3/bias/Regularizer/Const?
dense_3/bias/Regularizer/SumSum#dense_3/bias/Regularizer/Square:y:0'dense_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/Sum?
dense_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_3/bias/Regularizer/mul/x?
dense_3/bias/Regularizer/mulMul'dense_3/bias/Regularizer/mul/x:output:0%dense_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/mul?
IdentityIdentity(dense_3/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????22

Identity?
NoOpNoOp.^batch_normalization_3/StatefulPartitionedCall.^batch_normalization_4/StatefulPartitionedCall.^batch_normalization_5/StatefulPartitionedCall!^conv2d_2/StatefulPartitionedCall0^conv2d_2/bias/Regularizer/Square/ReadVariableOp2^conv2d_2/kernel/Regularizer/Square/ReadVariableOp!^conv2d_3/StatefulPartitionedCall0^conv2d_3/bias/Regularizer/Square/ReadVariableOp2^conv2d_3/kernel/Regularizer/Square/ReadVariableOp ^dense_2/StatefulPartitionedCall/^dense_2/bias/Regularizer/Square/ReadVariableOp1^dense_2/kernel/Regularizer/Square/ReadVariableOp ^dense_3/StatefulPartitionedCall/^dense_3/bias/Regularizer/Square/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*N
_input_shapes=
;:?????????d: : : : : : : : : : : : : : : : : : : : 2^
-batch_normalization_3/StatefulPartitionedCall-batch_normalization_3/StatefulPartitionedCall2^
-batch_normalization_4/StatefulPartitionedCall-batch_normalization_4/StatefulPartitionedCall2^
-batch_normalization_5/StatefulPartitionedCall-batch_normalization_5/StatefulPartitionedCall2D
 conv2d_2/StatefulPartitionedCall conv2d_2/StatefulPartitionedCall2b
/conv2d_2/bias/Regularizer/Square/ReadVariableOp/conv2d_2/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp1conv2d_2/kernel/Regularizer/Square/ReadVariableOp2D
 conv2d_3/StatefulPartitionedCall conv2d_3/StatefulPartitionedCall2b
/conv2d_3/bias/Regularizer/Square/ReadVariableOp/conv2d_3/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp1conv2d_3/kernel/Regularizer/Square/ReadVariableOp2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2`
.dense_2/bias/Regularizer/Square/ReadVariableOp.dense_2/bias/Regularizer/Square/ReadVariableOp2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2`
.dense_3/bias/Regularizer/Square/ReadVariableOp.dense_3/bias/Regularizer/Square/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp:V R
'
_output_shapes
:?????????d
'
_user_specified_namedense_2_input
?
?
7__inference_batch_normalization_3_layer_call_fn_5331235

inputs
unknown:	?
	unknown_0:	?
	unknown_1:	?
	unknown_2:	?
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_53292392
StatefulPartitionedCall|
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:??????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:??????????: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
K
/__inference_up_sampling2d_layer_call_fn_5331323

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 *J
_output_shapes8
6:4????????????????????????????????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *S
fNRL
J__inference_up_sampling2d_layer_call_and_return_conditional_losses_53293332
PartitionedCall?
IdentityIdentityPartitionedCall:output:0*
T0*J
_output_shapes8
6:4????????????????????????????????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4????????????????????????????????????:r n
J
_output_shapes8
6:4????????????????????????????????????
 
_user_specified_nameinputs
?	
?
7__inference_batch_normalization_5_layer_call_fn_5331637

inputs
unknown:	?
	unknown_0:	?
	unknown_1:	?
	unknown_2:	?
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_53298292
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*0
_output_shapes
:??????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:??????????: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:X T
0
_output_shapes
:??????????
 
_user_specified_nameinputs
?
d
F__inference_dropout_1_layer_call_and_return_conditional_losses_5329682

inputs

identity_1[
IdentityIdentityinputs*
T0*(
_output_shapes
:??????????2

Identityj

Identity_1IdentityIdentity:output:0*
T0*(
_output_shapes
:??????????2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:??????????:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5329764

inputs%
readvariableop_resource:@'
readvariableop_1_resource:@6
(fusedbatchnormv3_readvariableop_resource:@8
*fusedbatchnormv3_readvariableop_1_resource:@
identity??FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1t
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
:@*
dtype02
ReadVariableOpz
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes
:@*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes
:@*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes
:@*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*K
_output_shapes9
7:?????????@:@:@:@:@:*
epsilon%o?:*
is_training( 2
FusedBatchNormV3w
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*/
_output_shapes
:?????????@2

Identity?
NoOpNoOp ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*6
_input_shapes%
#:?????????@: : : : 2B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
?
b
F__inference_reshape_1_layer_call_and_return_conditional_losses_5329707

inputs
identityD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2?
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
strided_sliced
Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2
Reshape/shape/1d
Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :2
Reshape/shape/2d
Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :@2
Reshape/shape/3?
Reshape/shapePackstrided_slice:output:0Reshape/shape/1:output:0Reshape/shape/2:output:0Reshape/shape/3:output:0*
N*
T0*
_output_shapes
:2
Reshape/shapew
ReshapeReshapeinputsReshape/shape:output:0*
T0*/
_output_shapes
:?????????@2	
Reshapel
IdentityIdentityReshape:output:0*
T0*/
_output_shapes
:?????????@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:??????????:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
??
?
I__inference_sequential_1_layer_call_and_return_conditional_losses_5330292

inputs"
dense_2_5330191:	d?
dense_2_5330193:	?,
batch_normalization_3_5330197:	?,
batch_normalization_3_5330199:	?,
batch_normalization_3_5330201:	?,
batch_normalization_3_5330203:	?*
conv2d_2_5330207:@@
conv2d_2_5330209:@+
batch_normalization_4_5330213:@+
batch_normalization_4_5330215:@+
batch_normalization_4_5330217:@+
batch_normalization_4_5330219:@+
conv2d_3_5330222:@?
conv2d_3_5330224:	?,
batch_normalization_5_5330228:	?,
batch_normalization_5_5330230:	?,
batch_normalization_5_5330232:	?,
batch_normalization_5_5330234:	?"
dense_3_5330238:	?@2
dense_3_5330240:2
identity??-batch_normalization_3/StatefulPartitionedCall?-batch_normalization_4/StatefulPartitionedCall?-batch_normalization_5/StatefulPartitionedCall? conv2d_2/StatefulPartitionedCall?/conv2d_2/bias/Regularizer/Square/ReadVariableOp?1conv2d_2/kernel/Regularizer/Square/ReadVariableOp? conv2d_3/StatefulPartitionedCall?/conv2d_3/bias/Regularizer/Square/ReadVariableOp?1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?dense_2/StatefulPartitionedCall?.dense_2/bias/Regularizer/Square/ReadVariableOp?0dense_2/kernel/Regularizer/Square/ReadVariableOp?dense_3/StatefulPartitionedCall?.dense_3/bias/Regularizer/Square/ReadVariableOp?0dense_3/kernel/Regularizer/Square/ReadVariableOp?!dropout_1/StatefulPartitionedCall?
dense_2/StatefulPartitionedCallStatefulPartitionedCallinputsdense_2_5330191dense_2_5330193*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_2_layer_call_and_return_conditional_losses_53296712!
dense_2/StatefulPartitionedCall?
!dropout_1/StatefulPartitionedCallStatefulPartitionedCall(dense_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_dropout_1_layer_call_and_return_conditional_losses_53301282#
!dropout_1/StatefulPartitionedCall?
-batch_normalization_3/StatefulPartitionedCallStatefulPartitionedCall*dropout_1/StatefulPartitionedCall:output:0batch_normalization_3_5330197batch_normalization_3_5330199batch_normalization_3_5330201batch_normalization_3_5330203*
Tin	
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_53292392/
-batch_normalization_3/StatefulPartitionedCall?
reshape_1/PartitionedCallPartitionedCall6batch_normalization_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_reshape_1_layer_call_and_return_conditional_losses_53297072
reshape_1/PartitionedCall?
 conv2d_2/StatefulPartitionedCallStatefulPartitionedCall"reshape_1/PartitionedCall:output:0conv2d_2_5330207conv2d_2_5330209*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv2d_2_layer_call_and_return_conditional_losses_53297322"
 conv2d_2/StatefulPartitionedCall?
up_sampling2d/PartitionedCallPartitionedCall)conv2d_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *S
fNRL
J__inference_up_sampling2d_layer_call_and_return_conditional_losses_53297452
up_sampling2d/PartitionedCall?
-batch_normalization_4/StatefulPartitionedCallStatefulPartitionedCall&up_sampling2d/PartitionedCall:output:0batch_normalization_4_5330213batch_normalization_4_5330215batch_normalization_4_5330217batch_normalization_4_5330219*
Tin	
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_53300762/
-batch_normalization_4/StatefulPartitionedCall?
 conv2d_3/StatefulPartitionedCallStatefulPartitionedCall6batch_normalization_4/StatefulPartitionedCall:output:0conv2d_3_5330222conv2d_3_5330224*
Tin
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv2d_3_layer_call_and_return_conditional_losses_53297972"
 conv2d_3/StatefulPartitionedCall?
up_sampling2d_1/PartitionedCallPartitionedCall)conv2d_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_up_sampling2d_1_layer_call_and_return_conditional_losses_53298102!
up_sampling2d_1/PartitionedCall?
-batch_normalization_5/StatefulPartitionedCallStatefulPartitionedCall(up_sampling2d_1/PartitionedCall:output:0batch_normalization_5_5330228batch_normalization_5_5330230batch_normalization_5_5330232batch_normalization_5_5330234*
Tin	
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_53300172/
-batch_normalization_5/StatefulPartitionedCall?
flatten_1/PartitionedCallPartitionedCall6batch_normalization_5/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_flatten_1_layer_call_and_return_conditional_losses_53298452
flatten_1/PartitionedCall?
dense_3/StatefulPartitionedCallStatefulPartitionedCall"flatten_1/PartitionedCall:output:0dense_3_5330238dense_3_5330240*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????2*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_3_layer_call_and_return_conditional_losses_53298702!
dense_3/StatefulPartitionedCall?
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_5330191*
_output_shapes
:	d?*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp?
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	d?2#
!dense_2/kernel/Regularizer/Square?
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const?
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/Sum?
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_2/kernel/Regularizer/mul/x?
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mul?
.dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_5330193*
_output_shapes	
:?*
dtype020
.dense_2/bias/Regularizer/Square/ReadVariableOp?
dense_2/bias/Regularizer/SquareSquare6dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2!
dense_2/bias/Regularizer/Square?
dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_2/bias/Regularizer/Const?
dense_2/bias/Regularizer/SumSum#dense_2/bias/Regularizer/Square:y:0'dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/Sum?
dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_2/bias/Regularizer/mul/x?
dense_2/bias/Regularizer/mulMul'dense_2/bias/Regularizer/mul/x:output:0%dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/mul?
1conv2d_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_2_5330207*&
_output_shapes
:@@*
dtype023
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_2/kernel/Regularizer/SquareSquare9conv2d_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*&
_output_shapes
:@@2$
"conv2d_2/kernel/Regularizer/Square?
!conv2d_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_2/kernel/Regularizer/Const?
conv2d_2/kernel/Regularizer/SumSum&conv2d_2/kernel/Regularizer/Square:y:0*conv2d_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/Sum?
!conv2d_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_2/kernel/Regularizer/mul/x?
conv2d_2/kernel/Regularizer/mulMul*conv2d_2/kernel/Regularizer/mul/x:output:0(conv2d_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/mul?
/conv2d_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_2_5330209*
_output_shapes
:@*
dtype021
/conv2d_2/bias/Regularizer/Square/ReadVariableOp?
 conv2d_2/bias/Regularizer/SquareSquare7conv2d_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2"
 conv2d_2/bias/Regularizer/Square?
conv2d_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_2/bias/Regularizer/Const?
conv2d_2/bias/Regularizer/SumSum$conv2d_2/bias/Regularizer/Square:y:0(conv2d_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/Sum?
conv2d_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_2/bias/Regularizer/mul/x?
conv2d_2/bias/Regularizer/mulMul(conv2d_2/bias/Regularizer/mul/x:output:0&conv2d_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/mul?
1conv2d_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_3_5330222*'
_output_shapes
:@?*
dtype023
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_3/kernel/Regularizer/SquareSquare9conv2d_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*'
_output_shapes
:@?2$
"conv2d_3/kernel/Regularizer/Square?
!conv2d_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_3/kernel/Regularizer/Const?
conv2d_3/kernel/Regularizer/SumSum&conv2d_3/kernel/Regularizer/Square:y:0*conv2d_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/Sum?
!conv2d_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_3/kernel/Regularizer/mul/x?
conv2d_3/kernel/Regularizer/mulMul*conv2d_3/kernel/Regularizer/mul/x:output:0(conv2d_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/mul?
/conv2d_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_3_5330224*
_output_shapes	
:?*
dtype021
/conv2d_3/bias/Regularizer/Square/ReadVariableOp?
 conv2d_3/bias/Regularizer/SquareSquare7conv2d_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2"
 conv2d_3/bias/Regularizer/Square?
conv2d_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_3/bias/Regularizer/Const?
conv2d_3/bias/Regularizer/SumSum$conv2d_3/bias/Regularizer/Square:y:0(conv2d_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/Sum?
conv2d_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_3/bias/Regularizer/mul/x?
conv2d_3/bias/Regularizer/mulMul(conv2d_3/bias/Regularizer/mul/x:output:0&conv2d_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/mul?
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_5330238*
_output_shapes
:	?@2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp?
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?@22#
!dense_3/kernel/Regularizer/Square?
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const?
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum?
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_3/kernel/Regularizer/mul/x?
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul?
.dense_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_5330240*
_output_shapes
:2*
dtype020
.dense_3/bias/Regularizer/Square/ReadVariableOp?
dense_3/bias/Regularizer/SquareSquare6dense_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:22!
dense_3/bias/Regularizer/Square?
dense_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_3/bias/Regularizer/Const?
dense_3/bias/Regularizer/SumSum#dense_3/bias/Regularizer/Square:y:0'dense_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/Sum?
dense_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_3/bias/Regularizer/mul/x?
dense_3/bias/Regularizer/mulMul'dense_3/bias/Regularizer/mul/x:output:0%dense_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/mul?
IdentityIdentity(dense_3/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????22

Identity?
NoOpNoOp.^batch_normalization_3/StatefulPartitionedCall.^batch_normalization_4/StatefulPartitionedCall.^batch_normalization_5/StatefulPartitionedCall!^conv2d_2/StatefulPartitionedCall0^conv2d_2/bias/Regularizer/Square/ReadVariableOp2^conv2d_2/kernel/Regularizer/Square/ReadVariableOp!^conv2d_3/StatefulPartitionedCall0^conv2d_3/bias/Regularizer/Square/ReadVariableOp2^conv2d_3/kernel/Regularizer/Square/ReadVariableOp ^dense_2/StatefulPartitionedCall/^dense_2/bias/Regularizer/Square/ReadVariableOp1^dense_2/kernel/Regularizer/Square/ReadVariableOp ^dense_3/StatefulPartitionedCall/^dense_3/bias/Regularizer/Square/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp"^dropout_1/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*N
_input_shapes=
;:?????????d: : : : : : : : : : : : : : : : : : : : 2^
-batch_normalization_3/StatefulPartitionedCall-batch_normalization_3/StatefulPartitionedCall2^
-batch_normalization_4/StatefulPartitionedCall-batch_normalization_4/StatefulPartitionedCall2^
-batch_normalization_5/StatefulPartitionedCall-batch_normalization_5/StatefulPartitionedCall2D
 conv2d_2/StatefulPartitionedCall conv2d_2/StatefulPartitionedCall2b
/conv2d_2/bias/Regularizer/Square/ReadVariableOp/conv2d_2/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp1conv2d_2/kernel/Regularizer/Square/ReadVariableOp2D
 conv2d_3/StatefulPartitionedCall conv2d_3/StatefulPartitionedCall2b
/conv2d_3/bias/Regularizer/Square/ReadVariableOp/conv2d_3/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp1conv2d_3/kernel/Regularizer/Square/ReadVariableOp2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2`
.dense_2/bias/Regularizer/Square/ReadVariableOp.dense_2/bias/Regularizer/Square/ReadVariableOp2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2`
.dense_3/bias/Regularizer/Square/ReadVariableOp.dense_3/bias/Regularizer/Square/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2F
!dropout_1/StatefulPartitionedCall!dropout_1/StatefulPartitionedCall:O K
'
_output_shapes
:?????????d
 
_user_specified_nameinputs
?*
?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_5329239

inputs6
'assignmovingavg_readvariableop_resource:	?8
)assignmovingavg_1_readvariableop_resource:	?4
%batchnorm_mul_readvariableop_resource:	?0
!batchnorm_readvariableop_resource:	?
identity??AssignMovingAvg?AssignMovingAvg/ReadVariableOp?AssignMovingAvg_1? AssignMovingAvg_1/ReadVariableOp?batchnorm/ReadVariableOp?batchnorm/mul/ReadVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes
:	?*
	keep_dims(2
moments/mean}
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes
:	?2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*(
_output_shapes
:??????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes
:	?*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes	
:?*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes	
:?*
squeeze_dims
 2
moments/Squeeze_1s
AssignMovingAvg/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOp'assignmovingavg_readvariableop_resource*
_output_shapes	
:?*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes	
:?2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes	
:?2
AssignMovingAvg/mul?
AssignMovingAvgAssignSubVariableOp'assignmovingavg_readvariableop_resourceAssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp*
_output_shapes
 *
dtype02
AssignMovingAvgw
AssignMovingAvg_1/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOp)assignmovingavg_1_readvariableop_resource*
_output_shapes	
:?*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes	
:?2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes	
:?2
AssignMovingAvg_1/mul?
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype02
AssignMovingAvg_1g
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes	
:?2
batchnorm/addd
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes	
:?2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes	
:?*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2
batchnorm/mulw
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*(
_output_shapes
:??????????2
batchnorm/mul_1|
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes	
:?2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes	
:?*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes	
:?2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*(
_output_shapes
:??????????2
batchnorm/add_1o
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*(
_output_shapes
:??????????2

Identity?
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^batchnorm/ReadVariableOp^batchnorm/mul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:??????????: : : : 2"
AssignMovingAvgAssignMovingAvg2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp2&
AssignMovingAvg_1AssignMovingAvg_12D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp24
batchnorm/ReadVariableOpbatchnorm/ReadVariableOp2<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5329375

inputs%
readvariableop_resource:@'
readvariableop_1_resource:@6
(fusedbatchnormv3_readvariableop_resource:@8
*fusedbatchnormv3_readvariableop_1_resource:@
identity??FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1t
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
:@*
dtype02
ReadVariableOpz
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes
:@*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes
:@*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes
:@*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*]
_output_shapesK
I:+???????????????????????????@:@:@:@:@:*
epsilon%o?:*
is_training( 2
FusedBatchNormV3?
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*A
_output_shapes/
-:+???????????????????????????@2

Identity?
NoOpNoOp ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*H
_input_shapes7
5:+???????????????????????????@: : : : 2B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:i e
A
_output_shapes/
-:+???????????????????????????@
 
_user_specified_nameinputs
?
?
__inference_loss_fn_6_5331782L
9dense_3_kernel_regularizer_square_readvariableop_resource:	?@2
identity??0dense_3/kernel/Regularizer/Square/ReadVariableOp?
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp9dense_3_kernel_regularizer_square_readvariableop_resource*
_output_shapes
:	?@2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp?
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?@22#
!dense_3/kernel/Regularizer/Square?
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const?
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum?
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_3/kernel/Regularizer/mul/x?
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mull
IdentityIdentity"dense_3/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 2

Identity?
NoOpNoOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp
?
?
*__inference_conv2d_3_layer_call_fn_5331496

inputs"
unknown:@?
	unknown_0:	?
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv2d_3_layer_call_and_return_conditional_losses_53297972
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*0
_output_shapes
:??????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????@: : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
?Z
?
#__inference__traced_restore_5331946
file_prefix2
assignvariableop_dense_2_kernel:	d?.
assignvariableop_1_dense_2_bias:	?=
.assignvariableop_2_batch_normalization_3_gamma:	?<
-assignvariableop_3_batch_normalization_3_beta:	?C
4assignvariableop_4_batch_normalization_3_moving_mean:	?G
8assignvariableop_5_batch_normalization_3_moving_variance:	?<
"assignvariableop_6_conv2d_2_kernel:@@.
 assignvariableop_7_conv2d_2_bias:@<
.assignvariableop_8_batch_normalization_4_gamma:@;
-assignvariableop_9_batch_normalization_4_beta:@C
5assignvariableop_10_batch_normalization_4_moving_mean:@G
9assignvariableop_11_batch_normalization_4_moving_variance:@>
#assignvariableop_12_conv2d_3_kernel:@?0
!assignvariableop_13_conv2d_3_bias:	?>
/assignvariableop_14_batch_normalization_5_gamma:	?=
.assignvariableop_15_batch_normalization_5_beta:	?D
5assignvariableop_16_batch_normalization_5_moving_mean:	?H
9assignvariableop_17_batch_normalization_5_moving_variance:	?5
"assignvariableop_18_dense_3_kernel:	?@2.
 assignvariableop_19_dense_3_bias:2
identity_21??AssignVariableOp?AssignVariableOp_1?AssignVariableOp_10?AssignVariableOp_11?AssignVariableOp_12?AssignVariableOp_13?AssignVariableOp_14?AssignVariableOp_15?AssignVariableOp_16?AssignVariableOp_17?AssignVariableOp_18?AssignVariableOp_19?AssignVariableOp_2?AssignVariableOp_3?AssignVariableOp_4?AssignVariableOp_5?AssignVariableOp_6?AssignVariableOp_7?AssignVariableOp_8?AssignVariableOp_9?

RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:*
dtype0*?	
value?	B?	B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-1/gamma/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/beta/.ATTRIBUTES/VARIABLE_VALUEB;layer_with_weights-1/moving_mean/.ATTRIBUTES/VARIABLE_VALUEB?layer_with_weights-1/moving_variance/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-3/gamma/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/beta/.ATTRIBUTES/VARIABLE_VALUEB;layer_with_weights-3/moving_mean/.ATTRIBUTES/VARIABLE_VALUEB?layer_with_weights-3/moving_variance/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-4/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-4/bias/.ATTRIBUTES/VARIABLE_VALUEB5layer_with_weights-5/gamma/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-5/beta/.ATTRIBUTES/VARIABLE_VALUEB;layer_with_weights-5/moving_mean/.ATTRIBUTES/VARIABLE_VALUEB?layer_with_weights-5/moving_variance/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-6/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-6/bias/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
RestoreV2/tensor_names?
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:*
dtype0*=
value4B2B B B B B B B B B B B B B B B B B B B B B 2
RestoreV2/shape_and_slices?
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*h
_output_shapesV
T:::::::::::::::::::::*#
dtypes
22
	RestoreV2g
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:2

Identity?
AssignVariableOpAssignVariableOpassignvariableop_dense_2_kernelIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOpk

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:2

Identity_1?
AssignVariableOp_1AssignVariableOpassignvariableop_1_dense_2_biasIdentity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_1k

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:2

Identity_2?
AssignVariableOp_2AssignVariableOp.assignvariableop_2_batch_normalization_3_gammaIdentity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_2k

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:2

Identity_3?
AssignVariableOp_3AssignVariableOp-assignvariableop_3_batch_normalization_3_betaIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_3k

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:2

Identity_4?
AssignVariableOp_4AssignVariableOp4assignvariableop_4_batch_normalization_3_moving_meanIdentity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_4k

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:2

Identity_5?
AssignVariableOp_5AssignVariableOp8assignvariableop_5_batch_normalization_3_moving_varianceIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_5k

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:2

Identity_6?
AssignVariableOp_6AssignVariableOp"assignvariableop_6_conv2d_2_kernelIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_6k

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:2

Identity_7?
AssignVariableOp_7AssignVariableOp assignvariableop_7_conv2d_2_biasIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_7k

Identity_8IdentityRestoreV2:tensors:8"/device:CPU:0*
T0*
_output_shapes
:2

Identity_8?
AssignVariableOp_8AssignVariableOp.assignvariableop_8_batch_normalization_4_gammaIdentity_8:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_8k

Identity_9IdentityRestoreV2:tensors:9"/device:CPU:0*
T0*
_output_shapes
:2

Identity_9?
AssignVariableOp_9AssignVariableOp-assignvariableop_9_batch_normalization_4_betaIdentity_9:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_9n
Identity_10IdentityRestoreV2:tensors:10"/device:CPU:0*
T0*
_output_shapes
:2
Identity_10?
AssignVariableOp_10AssignVariableOp5assignvariableop_10_batch_normalization_4_moving_meanIdentity_10:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_10n
Identity_11IdentityRestoreV2:tensors:11"/device:CPU:0*
T0*
_output_shapes
:2
Identity_11?
AssignVariableOp_11AssignVariableOp9assignvariableop_11_batch_normalization_4_moving_varianceIdentity_11:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_11n
Identity_12IdentityRestoreV2:tensors:12"/device:CPU:0*
T0*
_output_shapes
:2
Identity_12?
AssignVariableOp_12AssignVariableOp#assignvariableop_12_conv2d_3_kernelIdentity_12:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_12n
Identity_13IdentityRestoreV2:tensors:13"/device:CPU:0*
T0*
_output_shapes
:2
Identity_13?
AssignVariableOp_13AssignVariableOp!assignvariableop_13_conv2d_3_biasIdentity_13:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_13n
Identity_14IdentityRestoreV2:tensors:14"/device:CPU:0*
T0*
_output_shapes
:2
Identity_14?
AssignVariableOp_14AssignVariableOp/assignvariableop_14_batch_normalization_5_gammaIdentity_14:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_14n
Identity_15IdentityRestoreV2:tensors:15"/device:CPU:0*
T0*
_output_shapes
:2
Identity_15?
AssignVariableOp_15AssignVariableOp.assignvariableop_15_batch_normalization_5_betaIdentity_15:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_15n
Identity_16IdentityRestoreV2:tensors:16"/device:CPU:0*
T0*
_output_shapes
:2
Identity_16?
AssignVariableOp_16AssignVariableOp5assignvariableop_16_batch_normalization_5_moving_meanIdentity_16:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_16n
Identity_17IdentityRestoreV2:tensors:17"/device:CPU:0*
T0*
_output_shapes
:2
Identity_17?
AssignVariableOp_17AssignVariableOp9assignvariableop_17_batch_normalization_5_moving_varianceIdentity_17:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_17n
Identity_18IdentityRestoreV2:tensors:18"/device:CPU:0*
T0*
_output_shapes
:2
Identity_18?
AssignVariableOp_18AssignVariableOp"assignvariableop_18_dense_3_kernelIdentity_18:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_18n
Identity_19IdentityRestoreV2:tensors:19"/device:CPU:0*
T0*
_output_shapes
:2
Identity_19?
AssignVariableOp_19AssignVariableOp assignvariableop_19_dense_3_biasIdentity_19:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_199
NoOpNoOp"/device:CPU:0*
_output_shapes
 2
NoOp?
Identity_20Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9^NoOp"/device:CPU:0*
T0*
_output_shapes
: 2
Identity_20f
Identity_21IdentityIdentity_20:output:0^NoOp_1*
T0*
_output_shapes
: 2
Identity_21?
NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_10^AssignVariableOp_11^AssignVariableOp_12^AssignVariableOp_13^AssignVariableOp_14^AssignVariableOp_15^AssignVariableOp_16^AssignVariableOp_17^AssignVariableOp_18^AssignVariableOp_19^AssignVariableOp_2^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^AssignVariableOp_8^AssignVariableOp_9*"
_acd_function_control_output(*
_output_shapes
 2
NoOp_1"#
identity_21Identity_21:output:0*=
_input_shapes,
*: : : : : : : : : : : : : : : : : : : : : 2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12*
AssignVariableOp_10AssignVariableOp_102*
AssignVariableOp_11AssignVariableOp_112*
AssignVariableOp_12AssignVariableOp_122*
AssignVariableOp_13AssignVariableOp_132*
AssignVariableOp_14AssignVariableOp_142*
AssignVariableOp_15AssignVariableOp_152*
AssignVariableOp_16AssignVariableOp_162*
AssignVariableOp_17AssignVariableOp_172*
AssignVariableOp_18AssignVariableOp_182*
AssignVariableOp_19AssignVariableOp_192(
AssignVariableOp_2AssignVariableOp_22(
AssignVariableOp_3AssignVariableOp_32(
AssignVariableOp_4AssignVariableOp_42(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_72(
AssignVariableOp_8AssignVariableOp_82(
AssignVariableOp_9AssignVariableOp_9:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
?
?
__inference_loss_fn_4_5331760U
:conv2d_3_kernel_regularizer_square_readvariableop_resource:@?
identity??1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?
1conv2d_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp:conv2d_3_kernel_regularizer_square_readvariableop_resource*'
_output_shapes
:@?*
dtype023
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_3/kernel/Regularizer/SquareSquare9conv2d_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*'
_output_shapes
:@?2$
"conv2d_3/kernel/Regularizer/Square?
!conv2d_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_3/kernel/Regularizer/Const?
conv2d_3/kernel/Regularizer/SumSum&conv2d_3/kernel/Regularizer/Square:y:0*conv2d_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/Sum?
!conv2d_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_3/kernel/Regularizer/mul/x?
conv2d_3/kernel/Regularizer/mulMul*conv2d_3/kernel/Regularizer/mul/x:output:0(conv2d_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/mulm
IdentityIdentity#conv2d_3/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 2

Identity?
NoOpNoOp2^conv2d_3/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2f
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp1conv2d_3/kernel/Regularizer/Square/ReadVariableOp
?
?
E__inference_conv2d_2_layer_call_and_return_conditional_losses_5329732

inputs8
conv2d_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity??BiasAdd/ReadVariableOp?Conv2D/ReadVariableOp?/conv2d_2/bias/Regularizer/Square/ReadVariableOp?1conv2d_2/kernel/Regularizer/Square/ReadVariableOp?
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:@@*
dtype02
Conv2D/ReadVariableOp?
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????@*
paddingSAME*
strides
2
Conv2D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????@2	
BiasAdd]
EluEluBiasAdd:output:0*
T0*/
_output_shapes
:?????????@2
Elu?
1conv2d_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:@@*
dtype023
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_2/kernel/Regularizer/SquareSquare9conv2d_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*&
_output_shapes
:@@2$
"conv2d_2/kernel/Regularizer/Square?
!conv2d_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_2/kernel/Regularizer/Const?
conv2d_2/kernel/Regularizer/SumSum&conv2d_2/kernel/Regularizer/Square:y:0*conv2d_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/Sum?
!conv2d_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_2/kernel/Regularizer/mul/x?
conv2d_2/kernel/Regularizer/mulMul*conv2d_2/kernel/Regularizer/mul/x:output:0(conv2d_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/mul?
/conv2d_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype021
/conv2d_2/bias/Regularizer/Square/ReadVariableOp?
 conv2d_2/bias/Regularizer/SquareSquare7conv2d_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2"
 conv2d_2/bias/Regularizer/Square?
conv2d_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_2/bias/Regularizer/Const?
conv2d_2/bias/Regularizer/SumSum$conv2d_2/bias/Regularizer/Square:y:0(conv2d_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/Sum?
conv2d_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_2/bias/Regularizer/mul/x?
conv2d_2/bias/Regularizer/mulMul(conv2d_2/bias/Regularizer/mul/x:output:0&conv2d_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/mult
IdentityIdentityElu:activations:0^NoOp*
T0*/
_output_shapes
:?????????@2

Identity?
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp0^conv2d_2/bias/Regularizer/Square/ReadVariableOp2^conv2d_2/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp2b
/conv2d_2/bias/Regularizer/Square/ReadVariableOp/conv2d_2/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp1conv2d_2/kernel/Regularizer/Square/ReadVariableOp:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
?
?
__inference_loss_fn_7_5331793E
7dense_3_bias_regularizer_square_readvariableop_resource:2
identity??.dense_3/bias/Regularizer/Square/ReadVariableOp?
.dense_3/bias/Regularizer/Square/ReadVariableOpReadVariableOp7dense_3_bias_regularizer_square_readvariableop_resource*
_output_shapes
:2*
dtype020
.dense_3/bias/Regularizer/Square/ReadVariableOp?
dense_3/bias/Regularizer/SquareSquare6dense_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:22!
dense_3/bias/Regularizer/Square?
dense_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_3/bias/Regularizer/Const?
dense_3/bias/Regularizer/SumSum#dense_3/bias/Regularizer/Square:y:0'dense_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/Sum?
dense_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_3/bias/Regularizer/mul/x?
dense_3/bias/Regularizer/mulMul'dense_3/bias/Regularizer/mul/x:output:0%dense_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/mulj
IdentityIdentity dense_3/bias/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 2

Identity
NoOpNoOp/^dense_3/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2`
.dense_3/bias/Regularizer/Square/ReadVariableOp.dense_3/bias/Regularizer/Square/ReadVariableOp
?
?
7__inference_batch_normalization_4_layer_call_fn_5331452

inputs
unknown:@
	unknown_0:@
	unknown_1:@
	unknown_2:@
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_53300762
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:?????????@2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*6
_input_shapes%
#:?????????@: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
?
?
__inference_loss_fn_1_5331727F
7dense_2_bias_regularizer_square_readvariableop_resource:	?
identity??.dense_2/bias/Regularizer/Square/ReadVariableOp?
.dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp7dense_2_bias_regularizer_square_readvariableop_resource*
_output_shapes	
:?*
dtype020
.dense_2/bias/Regularizer/Square/ReadVariableOp?
dense_2/bias/Regularizer/SquareSquare6dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2!
dense_2/bias/Regularizer/Square?
dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_2/bias/Regularizer/Const?
dense_2/bias/Regularizer/SumSum#dense_2/bias/Regularizer/Square:y:0'dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/Sum?
dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_2/bias/Regularizer/mul/x?
dense_2/bias/Regularizer/mulMul'dense_2/bias/Regularizer/mul/x:output:0%dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/mulj
IdentityIdentity dense_2/bias/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 2

Identity
NoOpNoOp/^dense_2/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2`
.dense_2/bias/Regularizer/Square/ReadVariableOp.dense_2/bias/Regularizer/Square/ReadVariableOp
?	
?
7__inference_batch_normalization_5_layer_call_fn_5331611

inputs
unknown:	?
	unknown_0:	?
	unknown_1:	?
	unknown_2:	?
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *B
_output_shapes0
.:,????????????????????????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_53295372
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*B
_output_shapes0
.:,????????????????????????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:,????????????????????????????: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:j f
B
_output_shapes0
.:,????????????????????????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5329537

inputs&
readvariableop_resource:	?(
readvariableop_1_resource:	?7
(fusedbatchnormv3_readvariableop_resource:	?9
*fusedbatchnormv3_readvariableop_1_resource:	?
identity??FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1u
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp{
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes	
:?*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes	
:?*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*b
_output_shapesP
N:,????????????????????????????:?:?:?:?:*
epsilon%o?:*
is_training( 2
FusedBatchNormV3?
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*B
_output_shapes0
.:,????????????????????????????2

Identity?
NoOpNoOp ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:,????????????????????????????: : : : 2B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:j f
B
_output_shapes0
.:,????????????????????????????
 
_user_specified_nameinputs
?
?
.__inference_sequential_1_layer_call_fn_5331084

inputs
unknown:	d?
	unknown_0:	?
	unknown_1:	?
	unknown_2:	?
	unknown_3:	?
	unknown_4:	?#
	unknown_5:@@
	unknown_6:@
	unknown_7:@
	unknown_8:@
	unknown_9:@

unknown_10:@%

unknown_11:@?

unknown_12:	?

unknown_13:	?

unknown_14:	?

unknown_15:	?

unknown_16:	?

unknown_17:	?@2

unknown_18:2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18* 
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????2*0
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_sequential_1_layer_call_and_return_conditional_losses_53302922
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????22

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*N
_input_shapes=
;:?????????d: : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:?????????d
 
_user_specified_nameinputs
?*
?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_5331209

inputs6
'assignmovingavg_readvariableop_resource:	?8
)assignmovingavg_1_readvariableop_resource:	?4
%batchnorm_mul_readvariableop_resource:	?0
!batchnorm_readvariableop_resource:	?
identity??AssignMovingAvg?AssignMovingAvg/ReadVariableOp?AssignMovingAvg_1? AssignMovingAvg_1/ReadVariableOp?batchnorm/ReadVariableOp?batchnorm/mul/ReadVariableOp?
moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 2 
moments/mean/reduction_indices?
moments/meanMeaninputs'moments/mean/reduction_indices:output:0*
T0*
_output_shapes
:	?*
	keep_dims(2
moments/mean}
moments/StopGradientStopGradientmoments/mean:output:0*
T0*
_output_shapes
:	?2
moments/StopGradient?
moments/SquaredDifferenceSquaredDifferenceinputsmoments/StopGradient:output:0*
T0*(
_output_shapes
:??????????2
moments/SquaredDifference?
"moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 2$
"moments/variance/reduction_indices?
moments/varianceMeanmoments/SquaredDifference:z:0+moments/variance/reduction_indices:output:0*
T0*
_output_shapes
:	?*
	keep_dims(2
moments/variance?
moments/SqueezeSqueezemoments/mean:output:0*
T0*
_output_shapes	
:?*
squeeze_dims
 2
moments/Squeeze?
moments/Squeeze_1Squeezemoments/variance:output:0*
T0*
_output_shapes	
:?*
squeeze_dims
 2
moments/Squeeze_1s
AssignMovingAvg/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg/decay?
AssignMovingAvg/ReadVariableOpReadVariableOp'assignmovingavg_readvariableop_resource*
_output_shapes	
:?*
dtype02 
AssignMovingAvg/ReadVariableOp?
AssignMovingAvg/subSub&AssignMovingAvg/ReadVariableOp:value:0moments/Squeeze:output:0*
T0*
_output_shapes	
:?2
AssignMovingAvg/sub?
AssignMovingAvg/mulMulAssignMovingAvg/sub:z:0AssignMovingAvg/decay:output:0*
T0*
_output_shapes	
:?2
AssignMovingAvg/mul?
AssignMovingAvgAssignSubVariableOp'assignmovingavg_readvariableop_resourceAssignMovingAvg/mul:z:0^AssignMovingAvg/ReadVariableOp*
_output_shapes
 *
dtype02
AssignMovingAvgw
AssignMovingAvg_1/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
?#<2
AssignMovingAvg_1/decay?
 AssignMovingAvg_1/ReadVariableOpReadVariableOp)assignmovingavg_1_readvariableop_resource*
_output_shapes	
:?*
dtype02"
 AssignMovingAvg_1/ReadVariableOp?
AssignMovingAvg_1/subSub(AssignMovingAvg_1/ReadVariableOp:value:0moments/Squeeze_1:output:0*
T0*
_output_shapes	
:?2
AssignMovingAvg_1/sub?
AssignMovingAvg_1/mulMulAssignMovingAvg_1/sub:z:0 AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes	
:?2
AssignMovingAvg_1/mul?
AssignMovingAvg_1AssignSubVariableOp)assignmovingavg_1_readvariableop_resourceAssignMovingAvg_1/mul:z:0!^AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype02
AssignMovingAvg_1g
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2moments/Squeeze_1:output:0batchnorm/add/y:output:0*
T0*
_output_shapes	
:?2
batchnorm/addd
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes	
:?2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes	
:?*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2
batchnorm/mulw
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*(
_output_shapes
:??????????2
batchnorm/mul_1|
batchnorm/mul_2Mulmoments/Squeeze:output:0batchnorm/mul:z:0*
T0*
_output_shapes	
:?2
batchnorm/mul_2?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes	
:?*
dtype02
batchnorm/ReadVariableOp?
batchnorm/subSub batchnorm/ReadVariableOp:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes	
:?2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*(
_output_shapes
:??????????2
batchnorm/add_1o
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*(
_output_shapes
:??????????2

Identity?
NoOpNoOp^AssignMovingAvg^AssignMovingAvg/ReadVariableOp^AssignMovingAvg_1!^AssignMovingAvg_1/ReadVariableOp^batchnorm/ReadVariableOp^batchnorm/mul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:??????????: : : : 2"
AssignMovingAvgAssignMovingAvg2@
AssignMovingAvg/ReadVariableOpAssignMovingAvg/ReadVariableOp2&
AssignMovingAvg_1AssignMovingAvg_12D
 AssignMovingAvg_1/ReadVariableOp AssignMovingAvg_1/ReadVariableOp24
batchnorm/ReadVariableOpbatchnorm/ReadVariableOp2<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
7__inference_batch_normalization_3_layer_call_fn_5331222

inputs
unknown:	?
	unknown_0:	?
	unknown_1:	?
	unknown_2:	?
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_53291792
StatefulPartitionedCall|
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*(
_output_shapes
:??????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:??????????: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
.__inference_sequential_1_layer_call_fn_5331039

inputs
unknown:	d?
	unknown_0:	?
	unknown_1:	?
	unknown_2:	?
	unknown_3:	?
	unknown_4:	?#
	unknown_5:@@
	unknown_6:@
	unknown_7:@
	unknown_8:@
	unknown_9:@

unknown_10:@%

unknown_11:@?

unknown_12:	?

unknown_13:	?

unknown_14:	?

unknown_15:	?

unknown_16:	?

unknown_17:	?@2

unknown_18:2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18* 
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????2*6
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_sequential_1_layer_call_and_return_conditional_losses_53299252
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????22

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*N
_input_shapes=
;:?????????d: : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:O K
'
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
?
D__inference_dense_2_layer_call_and_return_conditional_losses_5331119

inputs1
matmul_readvariableop_resource:	d?.
biasadd_readvariableop_resource:	?
identity??BiasAdd/ReadVariableOp?MatMul/ReadVariableOp?.dense_2/bias/Regularizer/Square/ReadVariableOp?0dense_2/kernel/Regularizer/Square/ReadVariableOp?
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	d?*
dtype02
MatMul/ReadVariableOpt
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2	
BiasAddV
EluEluBiasAdd:output:0*
T0*(
_output_shapes
:??????????2
Elu?
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	d?*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp?
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	d?2#
!dense_2/kernel/Regularizer/Square?
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const?
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/Sum?
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_2/kernel/Regularizer/mul/x?
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mul?
.dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes	
:?*
dtype020
.dense_2/bias/Regularizer/Square/ReadVariableOp?
dense_2/bias/Regularizer/SquareSquare6dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2!
dense_2/bias/Regularizer/Square?
dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_2/bias/Regularizer/Const?
dense_2/bias/Regularizer/SumSum#dense_2/bias/Regularizer/Square:y:0'dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/Sum?
dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_2/bias/Regularizer/mul/x?
dense_2/bias/Regularizer/mulMul'dense_2/bias/Regularizer/mul/x:output:0%dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/mulm
IdentityIdentityElu:activations:0^NoOp*
T0*(
_output_shapes
:??????????2

Identity?
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp/^dense_2/bias/Regularizer/Square/ReadVariableOp1^dense_2/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime**
_input_shapes
:?????????d: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp2`
.dense_2/bias/Regularizer/Square/ReadVariableOp.dense_2/bias/Regularizer/Square/ReadVariableOp2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5331562

inputs&
readvariableop_resource:	?(
readvariableop_1_resource:	?7
(fusedbatchnormv3_readvariableop_resource:	?9
*fusedbatchnormv3_readvariableop_1_resource:	?
identity??AssignNewValue?AssignNewValue_1?FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1u
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp{
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes	
:?*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes	
:?*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*b
_output_shapesP
N:,????????????????????????????:?:?:?:?:*
epsilon%o?:*
exponential_avg_factor%
?#<2
FusedBatchNormV3?
AssignNewValueAssignVariableOp(fusedbatchnormv3_readvariableop_resourceFusedBatchNormV3:batch_mean:0 ^FusedBatchNormV3/ReadVariableOp*
_output_shapes
 *
dtype02
AssignNewValue?
AssignNewValue_1AssignVariableOp*fusedbatchnormv3_readvariableop_1_resource!FusedBatchNormV3:batch_variance:0"^FusedBatchNormV3/ReadVariableOp_1*
_output_shapes
 *
dtype02
AssignNewValue_1?
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*B
_output_shapes0
.:,????????????????????????????2

Identity?
NoOpNoOp^AssignNewValue^AssignNewValue_1 ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:,????????????????????????????: : : : 2 
AssignNewValueAssignNewValue2$
AssignNewValue_1AssignNewValue_12B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:j f
B
_output_shapes0
.:,????????????????????????????
 
_user_specified_nameinputs
?
?
__inference_loss_fn_3_5331749F
8conv2d_2_bias_regularizer_square_readvariableop_resource:@
identity??/conv2d_2/bias/Regularizer/Square/ReadVariableOp?
/conv2d_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp8conv2d_2_bias_regularizer_square_readvariableop_resource*
_output_shapes
:@*
dtype021
/conv2d_2/bias/Regularizer/Square/ReadVariableOp?
 conv2d_2/bias/Regularizer/SquareSquare7conv2d_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2"
 conv2d_2/bias/Regularizer/Square?
conv2d_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_2/bias/Regularizer/Const?
conv2d_2/bias/Regularizer/SumSum$conv2d_2/bias/Regularizer/Square:y:0(conv2d_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/Sum?
conv2d_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_2/bias/Regularizer/mul/x?
conv2d_2/bias/Regularizer/mulMul(conv2d_2/bias/Regularizer/mul/x:output:0&conv2d_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/mulk
IdentityIdentity!conv2d_2/bias/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 2

Identity?
NoOpNoOp0^conv2d_2/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2b
/conv2d_2/bias/Regularizer/Square/ReadVariableOp/conv2d_2/bias/Regularizer/Square/ReadVariableOp
?
h
L__inference_up_sampling2d_1_layer_call_and_return_conditional_losses_5331508

inputs
identityD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2?
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2
strided_slice_
ConstConst*
_output_shapes
:*
dtype0*
valueB"      2
Const^
mulMulstrided_slice:output:0Const:output:0*
T0*
_output_shapes
:2
mul?
resize/ResizeNearestNeighborResizeNearestNeighborinputsmul:z:0*
T0*J
_output_shapes8
6:4????????????????????????????????????*
half_pixel_centers(2
resize/ResizeNearestNeighbor?
IdentityIdentity-resize/ResizeNearestNeighbor:resized_images:0*
T0*J
_output_shapes8
6:4????????????????????????????????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4????????????????????????????????????:r n
J
_output_shapes8
6:4????????????????????????????????????
 
_user_specified_nameinputs
?
?
.__inference_sequential_1_layer_call_fn_5330380
dense_2_input
unknown:	d?
	unknown_0:	?
	unknown_1:	?
	unknown_2:	?
	unknown_3:	?
	unknown_4:	?#
	unknown_5:@@
	unknown_6:@
	unknown_7:@
	unknown_8:@
	unknown_9:@

unknown_10:@%

unknown_11:@?

unknown_12:	?

unknown_13:	?

unknown_14:	?

unknown_15:	?

unknown_16:	?

unknown_17:	?@2

unknown_18:2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCalldense_2_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6	unknown_7	unknown_8	unknown_9
unknown_10
unknown_11
unknown_12
unknown_13
unknown_14
unknown_15
unknown_16
unknown_17
unknown_18* 
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????2*0
_read_only_resource_inputs
	
*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_sequential_1_layer_call_and_return_conditional_losses_53302922
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????22

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*N
_input_shapes=
;:?????????d: : : : : : : : : : : : : : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:V R
'
_output_shapes
:?????????d
'
_user_specified_namedense_2_input
?	
?
7__inference_batch_normalization_4_layer_call_fn_5331426

inputs
unknown:@
	unknown_0:@
	unknown_1:@
	unknown_2:@
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *A
_output_shapes/
-:+???????????????????????????@*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_53294192
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*A
_output_shapes/
-:+???????????????????????????@2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*H
_input_shapes7
5:+???????????????????????????@: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:i e
A
_output_shapes/
-:+???????????????????????????@
 
_user_specified_nameinputs
?
?
__inference_loss_fn_0_5331716L
9dense_2_kernel_regularizer_square_readvariableop_resource:	d?
identity??0dense_2/kernel/Regularizer/Square/ReadVariableOp?
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp9dense_2_kernel_regularizer_square_readvariableop_resource*
_output_shapes
:	d?*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp?
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	d?2#
!dense_2/kernel/Regularizer/Square?
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const?
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/Sum?
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_2/kernel/Regularizer/mul/x?
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mull
IdentityIdentity"dense_2/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 2

Identity?
NoOpNoOp1^dense_2/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp
?
?
)__inference_dense_3_layer_call_fn_5331705

inputs
unknown:	?@2
	unknown_0:2
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????2*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_3_layer_call_and_return_conditional_losses_53298702
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????22

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:??????????@: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:??????????@
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_5329179

inputs0
!batchnorm_readvariableop_resource:	?4
%batchnorm_mul_readvariableop_resource:	?2
#batchnorm_readvariableop_1_resource:	?2
#batchnorm_readvariableop_2_resource:	?
identity??batchnorm/ReadVariableOp?batchnorm/ReadVariableOp_1?batchnorm/ReadVariableOp_2?batchnorm/mul/ReadVariableOp?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes	
:?*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes	
:?2
batchnorm/addd
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes	
:?2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes	
:?*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2
batchnorm/mulw
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*(
_output_shapes
:??????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes	
:?*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes	
:?2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes	
:?*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes	
:?2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*(
_output_shapes
:??????????2
batchnorm/add_1o
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*(
_output_shapes
:??????????2

Identity?
NoOpNoOp^batchnorm/ReadVariableOp^batchnorm/ReadVariableOp_1^batchnorm/ReadVariableOp_2^batchnorm/mul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:??????????: : : : 24
batchnorm/ReadVariableOpbatchnorm/ReadVariableOp28
batchnorm/ReadVariableOp_1batchnorm/ReadVariableOp_128
batchnorm/ReadVariableOp_2batchnorm/ReadVariableOp_22<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
??
?
I__inference_sequential_1_layer_call_and_return_conditional_losses_5330588
dense_2_input"
dense_2_5330487:	d?
dense_2_5330489:	?,
batch_normalization_3_5330493:	?,
batch_normalization_3_5330495:	?,
batch_normalization_3_5330497:	?,
batch_normalization_3_5330499:	?*
conv2d_2_5330503:@@
conv2d_2_5330505:@+
batch_normalization_4_5330509:@+
batch_normalization_4_5330511:@+
batch_normalization_4_5330513:@+
batch_normalization_4_5330515:@+
conv2d_3_5330518:@?
conv2d_3_5330520:	?,
batch_normalization_5_5330524:	?,
batch_normalization_5_5330526:	?,
batch_normalization_5_5330528:	?,
batch_normalization_5_5330530:	?"
dense_3_5330534:	?@2
dense_3_5330536:2
identity??-batch_normalization_3/StatefulPartitionedCall?-batch_normalization_4/StatefulPartitionedCall?-batch_normalization_5/StatefulPartitionedCall? conv2d_2/StatefulPartitionedCall?/conv2d_2/bias/Regularizer/Square/ReadVariableOp?1conv2d_2/kernel/Regularizer/Square/ReadVariableOp? conv2d_3/StatefulPartitionedCall?/conv2d_3/bias/Regularizer/Square/ReadVariableOp?1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?dense_2/StatefulPartitionedCall?.dense_2/bias/Regularizer/Square/ReadVariableOp?0dense_2/kernel/Regularizer/Square/ReadVariableOp?dense_3/StatefulPartitionedCall?.dense_3/bias/Regularizer/Square/ReadVariableOp?0dense_3/kernel/Regularizer/Square/ReadVariableOp?!dropout_1/StatefulPartitionedCall?
dense_2/StatefulPartitionedCallStatefulPartitionedCalldense_2_inputdense_2_5330487dense_2_5330489*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_2_layer_call_and_return_conditional_losses_53296712!
dense_2/StatefulPartitionedCall?
!dropout_1/StatefulPartitionedCallStatefulPartitionedCall(dense_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_dropout_1_layer_call_and_return_conditional_losses_53301282#
!dropout_1/StatefulPartitionedCall?
-batch_normalization_3/StatefulPartitionedCallStatefulPartitionedCall*dropout_1/StatefulPartitionedCall:output:0batch_normalization_3_5330493batch_normalization_3_5330495batch_normalization_3_5330497batch_normalization_3_5330499*
Tin	
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_53292392/
-batch_normalization_3/StatefulPartitionedCall?
reshape_1/PartitionedCallPartitionedCall6batch_normalization_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_reshape_1_layer_call_and_return_conditional_losses_53297072
reshape_1/PartitionedCall?
 conv2d_2/StatefulPartitionedCallStatefulPartitionedCall"reshape_1/PartitionedCall:output:0conv2d_2_5330503conv2d_2_5330505*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv2d_2_layer_call_and_return_conditional_losses_53297322"
 conv2d_2/StatefulPartitionedCall?
up_sampling2d/PartitionedCallPartitionedCall)conv2d_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *S
fNRL
J__inference_up_sampling2d_layer_call_and_return_conditional_losses_53297452
up_sampling2d/PartitionedCall?
-batch_normalization_4/StatefulPartitionedCallStatefulPartitionedCall&up_sampling2d/PartitionedCall:output:0batch_normalization_4_5330509batch_normalization_4_5330511batch_normalization_4_5330513batch_normalization_4_5330515*
Tin	
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_53300762/
-batch_normalization_4/StatefulPartitionedCall?
 conv2d_3/StatefulPartitionedCallStatefulPartitionedCall6batch_normalization_4/StatefulPartitionedCall:output:0conv2d_3_5330518conv2d_3_5330520*
Tin
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv2d_3_layer_call_and_return_conditional_losses_53297972"
 conv2d_3/StatefulPartitionedCall?
up_sampling2d_1/PartitionedCallPartitionedCall)conv2d_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_up_sampling2d_1_layer_call_and_return_conditional_losses_53298102!
up_sampling2d_1/PartitionedCall?
-batch_normalization_5/StatefulPartitionedCallStatefulPartitionedCall(up_sampling2d_1/PartitionedCall:output:0batch_normalization_5_5330524batch_normalization_5_5330526batch_normalization_5_5330528batch_normalization_5_5330530*
Tin	
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_53300172/
-batch_normalization_5/StatefulPartitionedCall?
flatten_1/PartitionedCallPartitionedCall6batch_normalization_5/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_flatten_1_layer_call_and_return_conditional_losses_53298452
flatten_1/PartitionedCall?
dense_3/StatefulPartitionedCallStatefulPartitionedCall"flatten_1/PartitionedCall:output:0dense_3_5330534dense_3_5330536*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????2*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_3_layer_call_and_return_conditional_losses_53298702!
dense_3/StatefulPartitionedCall?
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_5330487*
_output_shapes
:	d?*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp?
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	d?2#
!dense_2/kernel/Regularizer/Square?
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const?
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/Sum?
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_2/kernel/Regularizer/mul/x?
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mul?
.dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_5330489*
_output_shapes	
:?*
dtype020
.dense_2/bias/Regularizer/Square/ReadVariableOp?
dense_2/bias/Regularizer/SquareSquare6dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2!
dense_2/bias/Regularizer/Square?
dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_2/bias/Regularizer/Const?
dense_2/bias/Regularizer/SumSum#dense_2/bias/Regularizer/Square:y:0'dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/Sum?
dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_2/bias/Regularizer/mul/x?
dense_2/bias/Regularizer/mulMul'dense_2/bias/Regularizer/mul/x:output:0%dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/mul?
1conv2d_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_2_5330503*&
_output_shapes
:@@*
dtype023
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_2/kernel/Regularizer/SquareSquare9conv2d_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*&
_output_shapes
:@@2$
"conv2d_2/kernel/Regularizer/Square?
!conv2d_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_2/kernel/Regularizer/Const?
conv2d_2/kernel/Regularizer/SumSum&conv2d_2/kernel/Regularizer/Square:y:0*conv2d_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/Sum?
!conv2d_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_2/kernel/Regularizer/mul/x?
conv2d_2/kernel/Regularizer/mulMul*conv2d_2/kernel/Regularizer/mul/x:output:0(conv2d_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/mul?
/conv2d_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_2_5330505*
_output_shapes
:@*
dtype021
/conv2d_2/bias/Regularizer/Square/ReadVariableOp?
 conv2d_2/bias/Regularizer/SquareSquare7conv2d_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2"
 conv2d_2/bias/Regularizer/Square?
conv2d_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_2/bias/Regularizer/Const?
conv2d_2/bias/Regularizer/SumSum$conv2d_2/bias/Regularizer/Square:y:0(conv2d_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/Sum?
conv2d_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_2/bias/Regularizer/mul/x?
conv2d_2/bias/Regularizer/mulMul(conv2d_2/bias/Regularizer/mul/x:output:0&conv2d_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/mul?
1conv2d_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_3_5330518*'
_output_shapes
:@?*
dtype023
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_3/kernel/Regularizer/SquareSquare9conv2d_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*'
_output_shapes
:@?2$
"conv2d_3/kernel/Regularizer/Square?
!conv2d_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_3/kernel/Regularizer/Const?
conv2d_3/kernel/Regularizer/SumSum&conv2d_3/kernel/Regularizer/Square:y:0*conv2d_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/Sum?
!conv2d_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_3/kernel/Regularizer/mul/x?
conv2d_3/kernel/Regularizer/mulMul*conv2d_3/kernel/Regularizer/mul/x:output:0(conv2d_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/mul?
/conv2d_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_3_5330520*
_output_shapes	
:?*
dtype021
/conv2d_3/bias/Regularizer/Square/ReadVariableOp?
 conv2d_3/bias/Regularizer/SquareSquare7conv2d_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2"
 conv2d_3/bias/Regularizer/Square?
conv2d_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_3/bias/Regularizer/Const?
conv2d_3/bias/Regularizer/SumSum$conv2d_3/bias/Regularizer/Square:y:0(conv2d_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/Sum?
conv2d_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_3/bias/Regularizer/mul/x?
conv2d_3/bias/Regularizer/mulMul(conv2d_3/bias/Regularizer/mul/x:output:0&conv2d_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/mul?
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_5330534*
_output_shapes
:	?@2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp?
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?@22#
!dense_3/kernel/Regularizer/Square?
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const?
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum?
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_3/kernel/Regularizer/mul/x?
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul?
.dense_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_5330536*
_output_shapes
:2*
dtype020
.dense_3/bias/Regularizer/Square/ReadVariableOp?
dense_3/bias/Regularizer/SquareSquare6dense_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:22!
dense_3/bias/Regularizer/Square?
dense_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_3/bias/Regularizer/Const?
dense_3/bias/Regularizer/SumSum#dense_3/bias/Regularizer/Square:y:0'dense_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/Sum?
dense_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_3/bias/Regularizer/mul/x?
dense_3/bias/Regularizer/mulMul'dense_3/bias/Regularizer/mul/x:output:0%dense_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/mul?
IdentityIdentity(dense_3/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????22

Identity?
NoOpNoOp.^batch_normalization_3/StatefulPartitionedCall.^batch_normalization_4/StatefulPartitionedCall.^batch_normalization_5/StatefulPartitionedCall!^conv2d_2/StatefulPartitionedCall0^conv2d_2/bias/Regularizer/Square/ReadVariableOp2^conv2d_2/kernel/Regularizer/Square/ReadVariableOp!^conv2d_3/StatefulPartitionedCall0^conv2d_3/bias/Regularizer/Square/ReadVariableOp2^conv2d_3/kernel/Regularizer/Square/ReadVariableOp ^dense_2/StatefulPartitionedCall/^dense_2/bias/Regularizer/Square/ReadVariableOp1^dense_2/kernel/Regularizer/Square/ReadVariableOp ^dense_3/StatefulPartitionedCall/^dense_3/bias/Regularizer/Square/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp"^dropout_1/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*N
_input_shapes=
;:?????????d: : : : : : : : : : : : : : : : : : : : 2^
-batch_normalization_3/StatefulPartitionedCall-batch_normalization_3/StatefulPartitionedCall2^
-batch_normalization_4/StatefulPartitionedCall-batch_normalization_4/StatefulPartitionedCall2^
-batch_normalization_5/StatefulPartitionedCall-batch_normalization_5/StatefulPartitionedCall2D
 conv2d_2/StatefulPartitionedCall conv2d_2/StatefulPartitionedCall2b
/conv2d_2/bias/Regularizer/Square/ReadVariableOp/conv2d_2/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp1conv2d_2/kernel/Regularizer/Square/ReadVariableOp2D
 conv2d_3/StatefulPartitionedCall conv2d_3/StatefulPartitionedCall2b
/conv2d_3/bias/Regularizer/Square/ReadVariableOp/conv2d_3/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp1conv2d_3/kernel/Regularizer/Square/ReadVariableOp2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2`
.dense_2/bias/Regularizer/Square/ReadVariableOp.dense_2/bias/Regularizer/Square/ReadVariableOp2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2`
.dense_3/bias/Regularizer/Square/ReadVariableOp.dense_3/bias/Regularizer/Square/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp2F
!dropout_1/StatefulPartitionedCall!dropout_1/StatefulPartitionedCall:V R
'
_output_shapes
:?????????d
'
_user_specified_namedense_2_input
?
f
J__inference_up_sampling2d_layer_call_and_return_conditional_losses_5329745

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"      2
Constc
Const_1Const*
_output_shapes
:*
dtype0*
valueB"      2	
Const_1X
mulMulConst:output:0Const_1:output:0*
T0*
_output_shapes
:2
mul?
resize/ResizeNearestNeighborResizeNearestNeighborinputsmul:z:0*
T0*/
_output_shapes
:?????????@*
half_pixel_centers(2
resize/ResizeNearestNeighbor?
IdentityIdentity-resize/ResizeNearestNeighbor:resized_images:0*
T0*/
_output_shapes
:?????????@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????@:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
?
h
L__inference_up_sampling2d_1_layer_call_and_return_conditional_losses_5329810

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"      2
Constc
Const_1Const*
_output_shapes
:*
dtype0*
valueB"      2	
Const_1X
mulMulConst:output:0Const_1:output:0*
T0*
_output_shapes
:2
mul?
resize/ResizeNearestNeighborResizeNearestNeighborinputsmul:z:0*
T0*0
_output_shapes
:??????????*
half_pixel_centers(2
resize/ResizeNearestNeighbor?
IdentityIdentity-resize/ResizeNearestNeighbor:resized_images:0*
T0*0
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:??????????:X T
0
_output_shapes
:??????????
 
_user_specified_nameinputs
ڋ
?
I__inference_sequential_1_layer_call_and_return_conditional_losses_5329925

inputs"
dense_2_5329672:	d?
dense_2_5329674:	?,
batch_normalization_3_5329684:	?,
batch_normalization_3_5329686:	?,
batch_normalization_3_5329688:	?,
batch_normalization_3_5329690:	?*
conv2d_2_5329733:@@
conv2d_2_5329735:@+
batch_normalization_4_5329765:@+
batch_normalization_4_5329767:@+
batch_normalization_4_5329769:@+
batch_normalization_4_5329771:@+
conv2d_3_5329798:@?
conv2d_3_5329800:	?,
batch_normalization_5_5329830:	?,
batch_normalization_5_5329832:	?,
batch_normalization_5_5329834:	?,
batch_normalization_5_5329836:	?"
dense_3_5329871:	?@2
dense_3_5329873:2
identity??-batch_normalization_3/StatefulPartitionedCall?-batch_normalization_4/StatefulPartitionedCall?-batch_normalization_5/StatefulPartitionedCall? conv2d_2/StatefulPartitionedCall?/conv2d_2/bias/Regularizer/Square/ReadVariableOp?1conv2d_2/kernel/Regularizer/Square/ReadVariableOp? conv2d_3/StatefulPartitionedCall?/conv2d_3/bias/Regularizer/Square/ReadVariableOp?1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?dense_2/StatefulPartitionedCall?.dense_2/bias/Regularizer/Square/ReadVariableOp?0dense_2/kernel/Regularizer/Square/ReadVariableOp?dense_3/StatefulPartitionedCall?.dense_3/bias/Regularizer/Square/ReadVariableOp?0dense_3/kernel/Regularizer/Square/ReadVariableOp?
dense_2/StatefulPartitionedCallStatefulPartitionedCallinputsdense_2_5329672dense_2_5329674*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_2_layer_call_and_return_conditional_losses_53296712!
dense_2/StatefulPartitionedCall?
dropout_1/PartitionedCallPartitionedCall(dense_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_dropout_1_layer_call_and_return_conditional_losses_53296822
dropout_1/PartitionedCall?
-batch_normalization_3/StatefulPartitionedCallStatefulPartitionedCall"dropout_1/PartitionedCall:output:0batch_normalization_3_5329684batch_normalization_3_5329686batch_normalization_3_5329688batch_normalization_3_5329690*
Tin	
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_53291792/
-batch_normalization_3/StatefulPartitionedCall?
reshape_1/PartitionedCallPartitionedCall6batch_normalization_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_reshape_1_layer_call_and_return_conditional_losses_53297072
reshape_1/PartitionedCall?
 conv2d_2/StatefulPartitionedCallStatefulPartitionedCall"reshape_1/PartitionedCall:output:0conv2d_2_5329733conv2d_2_5329735*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv2d_2_layer_call_and_return_conditional_losses_53297322"
 conv2d_2/StatefulPartitionedCall?
up_sampling2d/PartitionedCallPartitionedCall)conv2d_2/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *S
fNRL
J__inference_up_sampling2d_layer_call_and_return_conditional_losses_53297452
up_sampling2d/PartitionedCall?
-batch_normalization_4/StatefulPartitionedCallStatefulPartitionedCall&up_sampling2d/PartitionedCall:output:0batch_normalization_4_5329765batch_normalization_4_5329767batch_normalization_4_5329769batch_normalization_4_5329771*
Tin	
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_53297642/
-batch_normalization_4/StatefulPartitionedCall?
 conv2d_3/StatefulPartitionedCallStatefulPartitionedCall6batch_normalization_4/StatefulPartitionedCall:output:0conv2d_3_5329798conv2d_3_5329800*
Tin
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *N
fIRG
E__inference_conv2d_3_layer_call_and_return_conditional_losses_53297972"
 conv2d_3/StatefulPartitionedCall?
up_sampling2d_1/PartitionedCallPartitionedCall)conv2d_3/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_up_sampling2d_1_layer_call_and_return_conditional_losses_53298102!
up_sampling2d_1/PartitionedCall?
-batch_normalization_5/StatefulPartitionedCallStatefulPartitionedCall(up_sampling2d_1/PartitionedCall:output:0batch_normalization_5_5329830batch_normalization_5_5329832batch_normalization_5_5329834batch_normalization_5_5329836*
Tin	
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????*&
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_53298292/
-batch_normalization_5/StatefulPartitionedCall?
flatten_1/PartitionedCallPartitionedCall6batch_normalization_5/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 *(
_output_shapes
:??????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_flatten_1_layer_call_and_return_conditional_losses_53298452
flatten_1/PartitionedCall?
dense_3/StatefulPartitionedCallStatefulPartitionedCall"flatten_1/PartitionedCall:output:0dense_3_5329871dense_3_5329873*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????2*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_3_layer_call_and_return_conditional_losses_53298702!
dense_3/StatefulPartitionedCall?
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_5329672*
_output_shapes
:	d?*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp?
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	d?2#
!dense_2/kernel/Regularizer/Square?
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const?
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/Sum?
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_2/kernel/Regularizer/mul/x?
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mul?
.dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_2_5329674*
_output_shapes	
:?*
dtype020
.dense_2/bias/Regularizer/Square/ReadVariableOp?
dense_2/bias/Regularizer/SquareSquare6dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2!
dense_2/bias/Regularizer/Square?
dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_2/bias/Regularizer/Const?
dense_2/bias/Regularizer/SumSum#dense_2/bias/Regularizer/Square:y:0'dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/Sum?
dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_2/bias/Regularizer/mul/x?
dense_2/bias/Regularizer/mulMul'dense_2/bias/Regularizer/mul/x:output:0%dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/mul?
1conv2d_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_2_5329733*&
_output_shapes
:@@*
dtype023
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_2/kernel/Regularizer/SquareSquare9conv2d_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*&
_output_shapes
:@@2$
"conv2d_2/kernel/Regularizer/Square?
!conv2d_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_2/kernel/Regularizer/Const?
conv2d_2/kernel/Regularizer/SumSum&conv2d_2/kernel/Regularizer/Square:y:0*conv2d_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/Sum?
!conv2d_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_2/kernel/Regularizer/mul/x?
conv2d_2/kernel/Regularizer/mulMul*conv2d_2/kernel/Regularizer/mul/x:output:0(conv2d_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/mul?
/conv2d_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_2_5329735*
_output_shapes
:@*
dtype021
/conv2d_2/bias/Regularizer/Square/ReadVariableOp?
 conv2d_2/bias/Regularizer/SquareSquare7conv2d_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2"
 conv2d_2/bias/Regularizer/Square?
conv2d_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_2/bias/Regularizer/Const?
conv2d_2/bias/Regularizer/SumSum$conv2d_2/bias/Regularizer/Square:y:0(conv2d_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/Sum?
conv2d_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_2/bias/Regularizer/mul/x?
conv2d_2/bias/Regularizer/mulMul(conv2d_2/bias/Regularizer/mul/x:output:0&conv2d_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/mul?
1conv2d_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_3_5329798*'
_output_shapes
:@?*
dtype023
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_3/kernel/Regularizer/SquareSquare9conv2d_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*'
_output_shapes
:@?2$
"conv2d_3/kernel/Regularizer/Square?
!conv2d_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_3/kernel/Regularizer/Const?
conv2d_3/kernel/Regularizer/SumSum&conv2d_3/kernel/Regularizer/Square:y:0*conv2d_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/Sum?
!conv2d_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_3/kernel/Regularizer/mul/x?
conv2d_3/kernel/Regularizer/mulMul*conv2d_3/kernel/Regularizer/mul/x:output:0(conv2d_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/mul?
/conv2d_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_3_5329800*
_output_shapes	
:?*
dtype021
/conv2d_3/bias/Regularizer/Square/ReadVariableOp?
 conv2d_3/bias/Regularizer/SquareSquare7conv2d_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2"
 conv2d_3/bias/Regularizer/Square?
conv2d_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_3/bias/Regularizer/Const?
conv2d_3/bias/Regularizer/SumSum$conv2d_3/bias/Regularizer/Square:y:0(conv2d_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/Sum?
conv2d_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_3/bias/Regularizer/mul/x?
conv2d_3/bias/Regularizer/mulMul(conv2d_3/bias/Regularizer/mul/x:output:0&conv2d_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/mul?
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_5329871*
_output_shapes
:	?@2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp?
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?@22#
!dense_3/kernel/Regularizer/Square?
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const?
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum?
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_3/kernel/Regularizer/mul/x?
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul?
.dense_3/bias/Regularizer/Square/ReadVariableOpReadVariableOpdense_3_5329873*
_output_shapes
:2*
dtype020
.dense_3/bias/Regularizer/Square/ReadVariableOp?
dense_3/bias/Regularizer/SquareSquare6dense_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:22!
dense_3/bias/Regularizer/Square?
dense_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_3/bias/Regularizer/Const?
dense_3/bias/Regularizer/SumSum#dense_3/bias/Regularizer/Square:y:0'dense_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/Sum?
dense_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_3/bias/Regularizer/mul/x?
dense_3/bias/Regularizer/mulMul'dense_3/bias/Regularizer/mul/x:output:0%dense_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/mul?
IdentityIdentity(dense_3/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????22

Identity?
NoOpNoOp.^batch_normalization_3/StatefulPartitionedCall.^batch_normalization_4/StatefulPartitionedCall.^batch_normalization_5/StatefulPartitionedCall!^conv2d_2/StatefulPartitionedCall0^conv2d_2/bias/Regularizer/Square/ReadVariableOp2^conv2d_2/kernel/Regularizer/Square/ReadVariableOp!^conv2d_3/StatefulPartitionedCall0^conv2d_3/bias/Regularizer/Square/ReadVariableOp2^conv2d_3/kernel/Regularizer/Square/ReadVariableOp ^dense_2/StatefulPartitionedCall/^dense_2/bias/Regularizer/Square/ReadVariableOp1^dense_2/kernel/Regularizer/Square/ReadVariableOp ^dense_3/StatefulPartitionedCall/^dense_3/bias/Regularizer/Square/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*N
_input_shapes=
;:?????????d: : : : : : : : : : : : : : : : : : : : 2^
-batch_normalization_3/StatefulPartitionedCall-batch_normalization_3/StatefulPartitionedCall2^
-batch_normalization_4/StatefulPartitionedCall-batch_normalization_4/StatefulPartitionedCall2^
-batch_normalization_5/StatefulPartitionedCall-batch_normalization_5/StatefulPartitionedCall2D
 conv2d_2/StatefulPartitionedCall conv2d_2/StatefulPartitionedCall2b
/conv2d_2/bias/Regularizer/Square/ReadVariableOp/conv2d_2/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp1conv2d_2/kernel/Regularizer/Square/ReadVariableOp2D
 conv2d_3/StatefulPartitionedCall conv2d_3/StatefulPartitionedCall2b
/conv2d_3/bias/Regularizer/Square/ReadVariableOp/conv2d_3/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp1conv2d_3/kernel/Regularizer/Square/ReadVariableOp2B
dense_2/StatefulPartitionedCalldense_2/StatefulPartitionedCall2`
.dense_2/bias/Regularizer/Square/ReadVariableOp.dense_2/bias/Regularizer/Square/ReadVariableOp2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp2B
dense_3/StatefulPartitionedCalldense_3/StatefulPartitionedCall2`
.dense_3/bias/Regularizer/Square/ReadVariableOp.dense_3/bias/Regularizer/Square/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:?????????d
 
_user_specified_nameinputs
?	
?
7__inference_batch_normalization_5_layer_call_fn_5331650

inputs
unknown:	?
	unknown_0:	?
	unknown_1:	?
	unknown_2:	?
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2*
Tin	
2*
Tout
2*
_collective_manager_ids
 *0
_output_shapes
:??????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *[
fVRT
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_53300172
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*0
_output_shapes
:??????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:??????????: : : : 22
StatefulPartitionedCallStatefulPartitionedCall:X T
0
_output_shapes
:??????????
 
_user_specified_nameinputs
??
?
I__inference_sequential_1_layer_call_and_return_conditional_losses_5330994

inputs9
&dense_2_matmul_readvariableop_resource:	d?6
'dense_2_biasadd_readvariableop_resource:	?L
=batch_normalization_3_assignmovingavg_readvariableop_resource:	?N
?batch_normalization_3_assignmovingavg_1_readvariableop_resource:	?J
;batch_normalization_3_batchnorm_mul_readvariableop_resource:	?F
7batch_normalization_3_batchnorm_readvariableop_resource:	?A
'conv2d_2_conv2d_readvariableop_resource:@@6
(conv2d_2_biasadd_readvariableop_resource:@;
-batch_normalization_4_readvariableop_resource:@=
/batch_normalization_4_readvariableop_1_resource:@L
>batch_normalization_4_fusedbatchnormv3_readvariableop_resource:@N
@batch_normalization_4_fusedbatchnormv3_readvariableop_1_resource:@B
'conv2d_3_conv2d_readvariableop_resource:@?7
(conv2d_3_biasadd_readvariableop_resource:	?<
-batch_normalization_5_readvariableop_resource:	?>
/batch_normalization_5_readvariableop_1_resource:	?M
>batch_normalization_5_fusedbatchnormv3_readvariableop_resource:	?O
@batch_normalization_5_fusedbatchnormv3_readvariableop_1_resource:	?9
&dense_3_matmul_readvariableop_resource:	?@25
'dense_3_biasadd_readvariableop_resource:2
identity??%batch_normalization_3/AssignMovingAvg?4batch_normalization_3/AssignMovingAvg/ReadVariableOp?'batch_normalization_3/AssignMovingAvg_1?6batch_normalization_3/AssignMovingAvg_1/ReadVariableOp?.batch_normalization_3/batchnorm/ReadVariableOp?2batch_normalization_3/batchnorm/mul/ReadVariableOp?$batch_normalization_4/AssignNewValue?&batch_normalization_4/AssignNewValue_1?5batch_normalization_4/FusedBatchNormV3/ReadVariableOp?7batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1?$batch_normalization_4/ReadVariableOp?&batch_normalization_4/ReadVariableOp_1?$batch_normalization_5/AssignNewValue?&batch_normalization_5/AssignNewValue_1?5batch_normalization_5/FusedBatchNormV3/ReadVariableOp?7batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1?$batch_normalization_5/ReadVariableOp?&batch_normalization_5/ReadVariableOp_1?conv2d_2/BiasAdd/ReadVariableOp?conv2d_2/Conv2D/ReadVariableOp?/conv2d_2/bias/Regularizer/Square/ReadVariableOp?1conv2d_2/kernel/Regularizer/Square/ReadVariableOp?conv2d_3/BiasAdd/ReadVariableOp?conv2d_3/Conv2D/ReadVariableOp?/conv2d_3/bias/Regularizer/Square/ReadVariableOp?1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?dense_2/BiasAdd/ReadVariableOp?dense_2/MatMul/ReadVariableOp?.dense_2/bias/Regularizer/Square/ReadVariableOp?0dense_2/kernel/Regularizer/Square/ReadVariableOp?dense_3/BiasAdd/ReadVariableOp?dense_3/MatMul/ReadVariableOp?.dense_3/bias/Regularizer/Square/ReadVariableOp?0dense_3/kernel/Regularizer/Square/ReadVariableOp?
dense_2/MatMul/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes
:	d?*
dtype02
dense_2/MatMul/ReadVariableOp?
dense_2/MatMulMatMulinputs%dense_2/MatMul/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
dense_2/MatMul?
dense_2/BiasAdd/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02 
dense_2/BiasAdd/ReadVariableOp?
dense_2/BiasAddBiasAdddense_2/MatMul:product:0&dense_2/BiasAdd/ReadVariableOp:value:0*
T0*(
_output_shapes
:??????????2
dense_2/BiasAddn
dense_2/EluEludense_2/BiasAdd:output:0*
T0*(
_output_shapes
:??????????2
dense_2/Eluw
dropout_1/dropout/ConstConst*
_output_shapes
: *
dtype0*
valueB
 *n۶?2
dropout_1/dropout/Const?
dropout_1/dropout/MulMuldense_2/Elu:activations:0 dropout_1/dropout/Const:output:0*
T0*(
_output_shapes
:??????????2
dropout_1/dropout/Mul{
dropout_1/dropout/ShapeShapedense_2/Elu:activations:0*
T0*
_output_shapes
:2
dropout_1/dropout/Shape?
.dropout_1/dropout/random_uniform/RandomUniformRandomUniform dropout_1/dropout/Shape:output:0*
T0*(
_output_shapes
:??????????*
dtype020
.dropout_1/dropout/random_uniform/RandomUniform?
 dropout_1/dropout/GreaterEqual/yConst*
_output_shapes
: *
dtype0*
valueB
 *???>2"
 dropout_1/dropout/GreaterEqual/y?
dropout_1/dropout/GreaterEqualGreaterEqual7dropout_1/dropout/random_uniform/RandomUniform:output:0)dropout_1/dropout/GreaterEqual/y:output:0*
T0*(
_output_shapes
:??????????2 
dropout_1/dropout/GreaterEqual?
dropout_1/dropout/CastCast"dropout_1/dropout/GreaterEqual:z:0*

DstT0*

SrcT0
*(
_output_shapes
:??????????2
dropout_1/dropout/Cast?
dropout_1/dropout/Mul_1Muldropout_1/dropout/Mul:z:0dropout_1/dropout/Cast:y:0*
T0*(
_output_shapes
:??????????2
dropout_1/dropout/Mul_1?
4batch_normalization_3/moments/mean/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 26
4batch_normalization_3/moments/mean/reduction_indices?
"batch_normalization_3/moments/meanMeandropout_1/dropout/Mul_1:z:0=batch_normalization_3/moments/mean/reduction_indices:output:0*
T0*
_output_shapes
:	?*
	keep_dims(2$
"batch_normalization_3/moments/mean?
*batch_normalization_3/moments/StopGradientStopGradient+batch_normalization_3/moments/mean:output:0*
T0*
_output_shapes
:	?2,
*batch_normalization_3/moments/StopGradient?
/batch_normalization_3/moments/SquaredDifferenceSquaredDifferencedropout_1/dropout/Mul_1:z:03batch_normalization_3/moments/StopGradient:output:0*
T0*(
_output_shapes
:??????????21
/batch_normalization_3/moments/SquaredDifference?
8batch_normalization_3/moments/variance/reduction_indicesConst*
_output_shapes
:*
dtype0*
valueB: 2:
8batch_normalization_3/moments/variance/reduction_indices?
&batch_normalization_3/moments/varianceMean3batch_normalization_3/moments/SquaredDifference:z:0Abatch_normalization_3/moments/variance/reduction_indices:output:0*
T0*
_output_shapes
:	?*
	keep_dims(2(
&batch_normalization_3/moments/variance?
%batch_normalization_3/moments/SqueezeSqueeze+batch_normalization_3/moments/mean:output:0*
T0*
_output_shapes	
:?*
squeeze_dims
 2'
%batch_normalization_3/moments/Squeeze?
'batch_normalization_3/moments/Squeeze_1Squeeze/batch_normalization_3/moments/variance:output:0*
T0*
_output_shapes	
:?*
squeeze_dims
 2)
'batch_normalization_3/moments/Squeeze_1?
+batch_normalization_3/AssignMovingAvg/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
?#<2-
+batch_normalization_3/AssignMovingAvg/decay?
4batch_normalization_3/AssignMovingAvg/ReadVariableOpReadVariableOp=batch_normalization_3_assignmovingavg_readvariableop_resource*
_output_shapes	
:?*
dtype026
4batch_normalization_3/AssignMovingAvg/ReadVariableOp?
)batch_normalization_3/AssignMovingAvg/subSub<batch_normalization_3/AssignMovingAvg/ReadVariableOp:value:0.batch_normalization_3/moments/Squeeze:output:0*
T0*
_output_shapes	
:?2+
)batch_normalization_3/AssignMovingAvg/sub?
)batch_normalization_3/AssignMovingAvg/mulMul-batch_normalization_3/AssignMovingAvg/sub:z:04batch_normalization_3/AssignMovingAvg/decay:output:0*
T0*
_output_shapes	
:?2+
)batch_normalization_3/AssignMovingAvg/mul?
%batch_normalization_3/AssignMovingAvgAssignSubVariableOp=batch_normalization_3_assignmovingavg_readvariableop_resource-batch_normalization_3/AssignMovingAvg/mul:z:05^batch_normalization_3/AssignMovingAvg/ReadVariableOp*
_output_shapes
 *
dtype02'
%batch_normalization_3/AssignMovingAvg?
-batch_normalization_3/AssignMovingAvg_1/decayConst*
_output_shapes
: *
dtype0*
valueB
 *
?#<2/
-batch_normalization_3/AssignMovingAvg_1/decay?
6batch_normalization_3/AssignMovingAvg_1/ReadVariableOpReadVariableOp?batch_normalization_3_assignmovingavg_1_readvariableop_resource*
_output_shapes	
:?*
dtype028
6batch_normalization_3/AssignMovingAvg_1/ReadVariableOp?
+batch_normalization_3/AssignMovingAvg_1/subSub>batch_normalization_3/AssignMovingAvg_1/ReadVariableOp:value:00batch_normalization_3/moments/Squeeze_1:output:0*
T0*
_output_shapes	
:?2-
+batch_normalization_3/AssignMovingAvg_1/sub?
+batch_normalization_3/AssignMovingAvg_1/mulMul/batch_normalization_3/AssignMovingAvg_1/sub:z:06batch_normalization_3/AssignMovingAvg_1/decay:output:0*
T0*
_output_shapes	
:?2-
+batch_normalization_3/AssignMovingAvg_1/mul?
'batch_normalization_3/AssignMovingAvg_1AssignSubVariableOp?batch_normalization_3_assignmovingavg_1_readvariableop_resource/batch_normalization_3/AssignMovingAvg_1/mul:z:07^batch_normalization_3/AssignMovingAvg_1/ReadVariableOp*
_output_shapes
 *
dtype02)
'batch_normalization_3/AssignMovingAvg_1?
%batch_normalization_3/batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2'
%batch_normalization_3/batchnorm/add/y?
#batch_normalization_3/batchnorm/addAddV20batch_normalization_3/moments/Squeeze_1:output:0.batch_normalization_3/batchnorm/add/y:output:0*
T0*
_output_shapes	
:?2%
#batch_normalization_3/batchnorm/add?
%batch_normalization_3/batchnorm/RsqrtRsqrt'batch_normalization_3/batchnorm/add:z:0*
T0*
_output_shapes	
:?2'
%batch_normalization_3/batchnorm/Rsqrt?
2batch_normalization_3/batchnorm/mul/ReadVariableOpReadVariableOp;batch_normalization_3_batchnorm_mul_readvariableop_resource*
_output_shapes	
:?*
dtype024
2batch_normalization_3/batchnorm/mul/ReadVariableOp?
#batch_normalization_3/batchnorm/mulMul)batch_normalization_3/batchnorm/Rsqrt:y:0:batch_normalization_3/batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2%
#batch_normalization_3/batchnorm/mul?
%batch_normalization_3/batchnorm/mul_1Muldropout_1/dropout/Mul_1:z:0'batch_normalization_3/batchnorm/mul:z:0*
T0*(
_output_shapes
:??????????2'
%batch_normalization_3/batchnorm/mul_1?
%batch_normalization_3/batchnorm/mul_2Mul.batch_normalization_3/moments/Squeeze:output:0'batch_normalization_3/batchnorm/mul:z:0*
T0*
_output_shapes	
:?2'
%batch_normalization_3/batchnorm/mul_2?
.batch_normalization_3/batchnorm/ReadVariableOpReadVariableOp7batch_normalization_3_batchnorm_readvariableop_resource*
_output_shapes	
:?*
dtype020
.batch_normalization_3/batchnorm/ReadVariableOp?
#batch_normalization_3/batchnorm/subSub6batch_normalization_3/batchnorm/ReadVariableOp:value:0)batch_normalization_3/batchnorm/mul_2:z:0*
T0*
_output_shapes	
:?2%
#batch_normalization_3/batchnorm/sub?
%batch_normalization_3/batchnorm/add_1AddV2)batch_normalization_3/batchnorm/mul_1:z:0'batch_normalization_3/batchnorm/sub:z:0*
T0*(
_output_shapes
:??????????2'
%batch_normalization_3/batchnorm/add_1{
reshape_1/ShapeShape)batch_normalization_3/batchnorm/add_1:z:0*
T0*
_output_shapes
:2
reshape_1/Shape?
reshape_1/strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB: 2
reshape_1/strided_slice/stack?
reshape_1/strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2!
reshape_1/strided_slice/stack_1?
reshape_1/strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2!
reshape_1/strided_slice/stack_2?
reshape_1/strided_sliceStridedSlicereshape_1/Shape:output:0&reshape_1/strided_slice/stack:output:0(reshape_1/strided_slice/stack_1:output:0(reshape_1/strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
: *
shrink_axis_mask2
reshape_1/strided_slicex
reshape_1/Reshape/shape/1Const*
_output_shapes
: *
dtype0*
value	B :2
reshape_1/Reshape/shape/1x
reshape_1/Reshape/shape/2Const*
_output_shapes
: *
dtype0*
value	B :2
reshape_1/Reshape/shape/2x
reshape_1/Reshape/shape/3Const*
_output_shapes
: *
dtype0*
value	B :@2
reshape_1/Reshape/shape/3?
reshape_1/Reshape/shapePack reshape_1/strided_slice:output:0"reshape_1/Reshape/shape/1:output:0"reshape_1/Reshape/shape/2:output:0"reshape_1/Reshape/shape/3:output:0*
N*
T0*
_output_shapes
:2
reshape_1/Reshape/shape?
reshape_1/ReshapeReshape)batch_normalization_3/batchnorm/add_1:z:0 reshape_1/Reshape/shape:output:0*
T0*/
_output_shapes
:?????????@2
reshape_1/Reshape?
conv2d_2/Conv2D/ReadVariableOpReadVariableOp'conv2d_2_conv2d_readvariableop_resource*&
_output_shapes
:@@*
dtype02 
conv2d_2/Conv2D/ReadVariableOp?
conv2d_2/Conv2DConv2Dreshape_1/Reshape:output:0&conv2d_2/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????@*
paddingSAME*
strides
2
conv2d_2/Conv2D?
conv2d_2/BiasAdd/ReadVariableOpReadVariableOp(conv2d_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype02!
conv2d_2/BiasAdd/ReadVariableOp?
conv2d_2/BiasAddBiasAddconv2d_2/Conv2D:output:0'conv2d_2/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????@2
conv2d_2/BiasAddx
conv2d_2/EluEluconv2d_2/BiasAdd:output:0*
T0*/
_output_shapes
:?????????@2
conv2d_2/Elu{
up_sampling2d/ConstConst*
_output_shapes
:*
dtype0*
valueB"      2
up_sampling2d/Const
up_sampling2d/Const_1Const*
_output_shapes
:*
dtype0*
valueB"      2
up_sampling2d/Const_1?
up_sampling2d/mulMulup_sampling2d/Const:output:0up_sampling2d/Const_1:output:0*
T0*
_output_shapes
:2
up_sampling2d/mul?
*up_sampling2d/resize/ResizeNearestNeighborResizeNearestNeighborconv2d_2/Elu:activations:0up_sampling2d/mul:z:0*
T0*/
_output_shapes
:?????????@*
half_pixel_centers(2,
*up_sampling2d/resize/ResizeNearestNeighbor?
$batch_normalization_4/ReadVariableOpReadVariableOp-batch_normalization_4_readvariableop_resource*
_output_shapes
:@*
dtype02&
$batch_normalization_4/ReadVariableOp?
&batch_normalization_4/ReadVariableOp_1ReadVariableOp/batch_normalization_4_readvariableop_1_resource*
_output_shapes
:@*
dtype02(
&batch_normalization_4/ReadVariableOp_1?
5batch_normalization_4/FusedBatchNormV3/ReadVariableOpReadVariableOp>batch_normalization_4_fusedbatchnormv3_readvariableop_resource*
_output_shapes
:@*
dtype027
5batch_normalization_4/FusedBatchNormV3/ReadVariableOp?
7batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1ReadVariableOp@batch_normalization_4_fusedbatchnormv3_readvariableop_1_resource*
_output_shapes
:@*
dtype029
7batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1?
&batch_normalization_4/FusedBatchNormV3FusedBatchNormV3;up_sampling2d/resize/ResizeNearestNeighbor:resized_images:0,batch_normalization_4/ReadVariableOp:value:0.batch_normalization_4/ReadVariableOp_1:value:0=batch_normalization_4/FusedBatchNormV3/ReadVariableOp:value:0?batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*K
_output_shapes9
7:?????????@:@:@:@:@:*
epsilon%o?:*
exponential_avg_factor%
?#<2(
&batch_normalization_4/FusedBatchNormV3?
$batch_normalization_4/AssignNewValueAssignVariableOp>batch_normalization_4_fusedbatchnormv3_readvariableop_resource3batch_normalization_4/FusedBatchNormV3:batch_mean:06^batch_normalization_4/FusedBatchNormV3/ReadVariableOp*
_output_shapes
 *
dtype02&
$batch_normalization_4/AssignNewValue?
&batch_normalization_4/AssignNewValue_1AssignVariableOp@batch_normalization_4_fusedbatchnormv3_readvariableop_1_resource7batch_normalization_4/FusedBatchNormV3:batch_variance:08^batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1*
_output_shapes
 *
dtype02(
&batch_normalization_4/AssignNewValue_1?
conv2d_3/Conv2D/ReadVariableOpReadVariableOp'conv2d_3_conv2d_readvariableop_resource*'
_output_shapes
:@?*
dtype02 
conv2d_3/Conv2D/ReadVariableOp?
conv2d_3/Conv2DConv2D*batch_normalization_4/FusedBatchNormV3:y:0&conv2d_3/Conv2D/ReadVariableOp:value:0*
T0*0
_output_shapes
:??????????*
paddingSAME*
strides
2
conv2d_3/Conv2D?
conv2d_3/BiasAdd/ReadVariableOpReadVariableOp(conv2d_3_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype02!
conv2d_3/BiasAdd/ReadVariableOp?
conv2d_3/BiasAddBiasAddconv2d_3/Conv2D:output:0'conv2d_3/BiasAdd/ReadVariableOp:value:0*
T0*0
_output_shapes
:??????????2
conv2d_3/BiasAddy
conv2d_3/EluEluconv2d_3/BiasAdd:output:0*
T0*0
_output_shapes
:??????????2
conv2d_3/Elu
up_sampling2d_1/ConstConst*
_output_shapes
:*
dtype0*
valueB"      2
up_sampling2d_1/Const?
up_sampling2d_1/Const_1Const*
_output_shapes
:*
dtype0*
valueB"      2
up_sampling2d_1/Const_1?
up_sampling2d_1/mulMulup_sampling2d_1/Const:output:0 up_sampling2d_1/Const_1:output:0*
T0*
_output_shapes
:2
up_sampling2d_1/mul?
,up_sampling2d_1/resize/ResizeNearestNeighborResizeNearestNeighborconv2d_3/Elu:activations:0up_sampling2d_1/mul:z:0*
T0*0
_output_shapes
:??????????*
half_pixel_centers(2.
,up_sampling2d_1/resize/ResizeNearestNeighbor?
$batch_normalization_5/ReadVariableOpReadVariableOp-batch_normalization_5_readvariableop_resource*
_output_shapes	
:?*
dtype02&
$batch_normalization_5/ReadVariableOp?
&batch_normalization_5/ReadVariableOp_1ReadVariableOp/batch_normalization_5_readvariableop_1_resource*
_output_shapes	
:?*
dtype02(
&batch_normalization_5/ReadVariableOp_1?
5batch_normalization_5/FusedBatchNormV3/ReadVariableOpReadVariableOp>batch_normalization_5_fusedbatchnormv3_readvariableop_resource*
_output_shapes	
:?*
dtype027
5batch_normalization_5/FusedBatchNormV3/ReadVariableOp?
7batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1ReadVariableOp@batch_normalization_5_fusedbatchnormv3_readvariableop_1_resource*
_output_shapes	
:?*
dtype029
7batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1?
&batch_normalization_5/FusedBatchNormV3FusedBatchNormV3=up_sampling2d_1/resize/ResizeNearestNeighbor:resized_images:0,batch_normalization_5/ReadVariableOp:value:0.batch_normalization_5/ReadVariableOp_1:value:0=batch_normalization_5/FusedBatchNormV3/ReadVariableOp:value:0?batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*P
_output_shapes>
<:??????????:?:?:?:?:*
epsilon%o?:*
exponential_avg_factor%
?#<2(
&batch_normalization_5/FusedBatchNormV3?
$batch_normalization_5/AssignNewValueAssignVariableOp>batch_normalization_5_fusedbatchnormv3_readvariableop_resource3batch_normalization_5/FusedBatchNormV3:batch_mean:06^batch_normalization_5/FusedBatchNormV3/ReadVariableOp*
_output_shapes
 *
dtype02&
$batch_normalization_5/AssignNewValue?
&batch_normalization_5/AssignNewValue_1AssignVariableOp@batch_normalization_5_fusedbatchnormv3_readvariableop_1_resource7batch_normalization_5/FusedBatchNormV3:batch_variance:08^batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1*
_output_shapes
 *
dtype02(
&batch_normalization_5/AssignNewValue_1s
flatten_1/ConstConst*
_output_shapes
:*
dtype0*
valueB"????    2
flatten_1/Const?
flatten_1/ReshapeReshape*batch_normalization_5/FusedBatchNormV3:y:0flatten_1/Const:output:0*
T0*(
_output_shapes
:??????????@2
flatten_1/Reshape?
dense_3/MatMul/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes
:	?@2*
dtype02
dense_3/MatMul/ReadVariableOp?
dense_3/MatMulMatMulflatten_1/Reshape:output:0%dense_3/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????22
dense_3/MatMul?
dense_3/BiasAdd/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype02 
dense_3/BiasAdd/ReadVariableOp?
dense_3/BiasAddBiasAdddense_3/MatMul:product:0&dense_3/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????22
dense_3/BiasAddy
dense_3/SigmoidSigmoiddense_3/BiasAdd:output:0*
T0*'
_output_shapes
:?????????22
dense_3/Sigmoid?
0dense_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_2_matmul_readvariableop_resource*
_output_shapes
:	d?*
dtype022
0dense_2/kernel/Regularizer/Square/ReadVariableOp?
!dense_2/kernel/Regularizer/SquareSquare8dense_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	d?2#
!dense_2/kernel/Regularizer/Square?
 dense_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_2/kernel/Regularizer/Const?
dense_2/kernel/Regularizer/SumSum%dense_2/kernel/Regularizer/Square:y:0)dense_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/Sum?
 dense_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_2/kernel/Regularizer/mul/x?
dense_2/kernel/Regularizer/mulMul)dense_2/kernel/Regularizer/mul/x:output:0'dense_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_2/kernel/Regularizer/mul?
.dense_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp'dense_2_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype020
.dense_2/bias/Regularizer/Square/ReadVariableOp?
dense_2/bias/Regularizer/SquareSquare6dense_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2!
dense_2/bias/Regularizer/Square?
dense_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_2/bias/Regularizer/Const?
dense_2/bias/Regularizer/SumSum#dense_2/bias/Regularizer/Square:y:0'dense_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/Sum?
dense_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_2/bias/Regularizer/mul/x?
dense_2/bias/Regularizer/mulMul'dense_2/bias/Regularizer/mul/x:output:0%dense_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_2/bias/Regularizer/mul?
1conv2d_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'conv2d_2_conv2d_readvariableop_resource*&
_output_shapes
:@@*
dtype023
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_2/kernel/Regularizer/SquareSquare9conv2d_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*&
_output_shapes
:@@2$
"conv2d_2/kernel/Regularizer/Square?
!conv2d_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_2/kernel/Regularizer/Const?
conv2d_2/kernel/Regularizer/SumSum&conv2d_2/kernel/Regularizer/Square:y:0*conv2d_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/Sum?
!conv2d_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_2/kernel/Regularizer/mul/x?
conv2d_2/kernel/Regularizer/mulMul*conv2d_2/kernel/Regularizer/mul/x:output:0(conv2d_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/mul?
/conv2d_2/bias/Regularizer/Square/ReadVariableOpReadVariableOp(conv2d_2_biasadd_readvariableop_resource*
_output_shapes
:@*
dtype021
/conv2d_2/bias/Regularizer/Square/ReadVariableOp?
 conv2d_2/bias/Regularizer/SquareSquare7conv2d_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2"
 conv2d_2/bias/Regularizer/Square?
conv2d_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_2/bias/Regularizer/Const?
conv2d_2/bias/Regularizer/SumSum$conv2d_2/bias/Regularizer/Square:y:0(conv2d_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/Sum?
conv2d_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_2/bias/Regularizer/mul/x?
conv2d_2/bias/Regularizer/mulMul(conv2d_2/bias/Regularizer/mul/x:output:0&conv2d_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/mul?
1conv2d_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp'conv2d_3_conv2d_readvariableop_resource*'
_output_shapes
:@?*
dtype023
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_3/kernel/Regularizer/SquareSquare9conv2d_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*'
_output_shapes
:@?2$
"conv2d_3/kernel/Regularizer/Square?
!conv2d_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_3/kernel/Regularizer/Const?
conv2d_3/kernel/Regularizer/SumSum&conv2d_3/kernel/Regularizer/Square:y:0*conv2d_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/Sum?
!conv2d_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_3/kernel/Regularizer/mul/x?
conv2d_3/kernel/Regularizer/mulMul*conv2d_3/kernel/Regularizer/mul/x:output:0(conv2d_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_3/kernel/Regularizer/mul?
/conv2d_3/bias/Regularizer/Square/ReadVariableOpReadVariableOp(conv2d_3_biasadd_readvariableop_resource*
_output_shapes	
:?*
dtype021
/conv2d_3/bias/Regularizer/Square/ReadVariableOp?
 conv2d_3/bias/Regularizer/SquareSquare7conv2d_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2"
 conv2d_3/bias/Regularizer/Square?
conv2d_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_3/bias/Regularizer/Const?
conv2d_3/bias/Regularizer/SumSum$conv2d_3/bias/Regularizer/Square:y:0(conv2d_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/Sum?
conv2d_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_3/bias/Regularizer/mul/x?
conv2d_3/bias/Regularizer/mulMul(conv2d_3/bias/Regularizer/mul/x:output:0&conv2d_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/mul?
0dense_3/kernel/Regularizer/Square/ReadVariableOpReadVariableOp&dense_3_matmul_readvariableop_resource*
_output_shapes
:	?@2*
dtype022
0dense_3/kernel/Regularizer/Square/ReadVariableOp?
!dense_3/kernel/Regularizer/SquareSquare8dense_3/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:	?@22#
!dense_3/kernel/Regularizer/Square?
 dense_3/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB"       2"
 dense_3/kernel/Regularizer/Const?
dense_3/kernel/Regularizer/SumSum%dense_3/kernel/Regularizer/Square:y:0)dense_3/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/Sum?
 dense_3/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2"
 dense_3/kernel/Regularizer/mul/x?
dense_3/kernel/Regularizer/mulMul)dense_3/kernel/Regularizer/mul/x:output:0'dense_3/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2 
dense_3/kernel/Regularizer/mul?
.dense_3/bias/Regularizer/Square/ReadVariableOpReadVariableOp'dense_3_biasadd_readvariableop_resource*
_output_shapes
:2*
dtype020
.dense_3/bias/Regularizer/Square/ReadVariableOp?
dense_3/bias/Regularizer/SquareSquare6dense_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:22!
dense_3/bias/Regularizer/Square?
dense_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2 
dense_3/bias/Regularizer/Const?
dense_3/bias/Regularizer/SumSum#dense_3/bias/Regularizer/Square:y:0'dense_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/Sum?
dense_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2 
dense_3/bias/Regularizer/mul/x?
dense_3/bias/Regularizer/mulMul'dense_3/bias/Regularizer/mul/x:output:0%dense_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
dense_3/bias/Regularizer/muln
IdentityIdentitydense_3/Sigmoid:y:0^NoOp*
T0*'
_output_shapes
:?????????22

Identity?
NoOpNoOp&^batch_normalization_3/AssignMovingAvg5^batch_normalization_3/AssignMovingAvg/ReadVariableOp(^batch_normalization_3/AssignMovingAvg_17^batch_normalization_3/AssignMovingAvg_1/ReadVariableOp/^batch_normalization_3/batchnorm/ReadVariableOp3^batch_normalization_3/batchnorm/mul/ReadVariableOp%^batch_normalization_4/AssignNewValue'^batch_normalization_4/AssignNewValue_16^batch_normalization_4/FusedBatchNormV3/ReadVariableOp8^batch_normalization_4/FusedBatchNormV3/ReadVariableOp_1%^batch_normalization_4/ReadVariableOp'^batch_normalization_4/ReadVariableOp_1%^batch_normalization_5/AssignNewValue'^batch_normalization_5/AssignNewValue_16^batch_normalization_5/FusedBatchNormV3/ReadVariableOp8^batch_normalization_5/FusedBatchNormV3/ReadVariableOp_1%^batch_normalization_5/ReadVariableOp'^batch_normalization_5/ReadVariableOp_1 ^conv2d_2/BiasAdd/ReadVariableOp^conv2d_2/Conv2D/ReadVariableOp0^conv2d_2/bias/Regularizer/Square/ReadVariableOp2^conv2d_2/kernel/Regularizer/Square/ReadVariableOp ^conv2d_3/BiasAdd/ReadVariableOp^conv2d_3/Conv2D/ReadVariableOp0^conv2d_3/bias/Regularizer/Square/ReadVariableOp2^conv2d_3/kernel/Regularizer/Square/ReadVariableOp^dense_2/BiasAdd/ReadVariableOp^dense_2/MatMul/ReadVariableOp/^dense_2/bias/Regularizer/Square/ReadVariableOp1^dense_2/kernel/Regularizer/Square/ReadVariableOp^dense_3/BiasAdd/ReadVariableOp^dense_3/MatMul/ReadVariableOp/^dense_3/bias/Regularizer/Square/ReadVariableOp1^dense_3/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*N
_input_shapes=
;:?????????d: : : : : : : : : : : : : : : : : : : : 2N
%batch_normalization_3/AssignMovingAvg%batch_normalization_3/AssignMovingAvg2l
4batch_normalization_3/AssignMovingAvg/ReadVariableOp4batch_normalization_3/AssignMovingAvg/ReadVariableOp2R
'batch_normalization_3/AssignMovingAvg_1'batch_normalization_3/AssignMovingAvg_12p
6batch_normalization_3/AssignMovingAvg_1/ReadVariableOp6batch_normalization_3/AssignMovingAvg_1/ReadVariableOp2`
.batch_normalization_3/batchnorm/ReadVariableOp.batch_normalization_3/batchnorm/ReadVariableOp2h
2batch_normalization_3/batchnorm/mul/ReadVariableOp2batch_normalization_3/batchnorm/mul/ReadVariableOp2L
$batch_normalization_4/AssignNewValue$batch_normalization_4/AssignNewValue2P
&batch_normalization_4/AssignNewValue_1&batch_normalization_4/AssignNewValue_12n
5batch_normalization_4/FusedBatchNormV3/ReadVariableOp5batch_normalization_4/FusedBatchNormV3/ReadVariableOp2r
7batch_normalization_4/FusedBatchNormV3/ReadVariableOp_17batch_normalization_4/FusedBatchNormV3/ReadVariableOp_12L
$batch_normalization_4/ReadVariableOp$batch_normalization_4/ReadVariableOp2P
&batch_normalization_4/ReadVariableOp_1&batch_normalization_4/ReadVariableOp_12L
$batch_normalization_5/AssignNewValue$batch_normalization_5/AssignNewValue2P
&batch_normalization_5/AssignNewValue_1&batch_normalization_5/AssignNewValue_12n
5batch_normalization_5/FusedBatchNormV3/ReadVariableOp5batch_normalization_5/FusedBatchNormV3/ReadVariableOp2r
7batch_normalization_5/FusedBatchNormV3/ReadVariableOp_17batch_normalization_5/FusedBatchNormV3/ReadVariableOp_12L
$batch_normalization_5/ReadVariableOp$batch_normalization_5/ReadVariableOp2P
&batch_normalization_5/ReadVariableOp_1&batch_normalization_5/ReadVariableOp_12B
conv2d_2/BiasAdd/ReadVariableOpconv2d_2/BiasAdd/ReadVariableOp2@
conv2d_2/Conv2D/ReadVariableOpconv2d_2/Conv2D/ReadVariableOp2b
/conv2d_2/bias/Regularizer/Square/ReadVariableOp/conv2d_2/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp1conv2d_2/kernel/Regularizer/Square/ReadVariableOp2B
conv2d_3/BiasAdd/ReadVariableOpconv2d_3/BiasAdd/ReadVariableOp2@
conv2d_3/Conv2D/ReadVariableOpconv2d_3/Conv2D/ReadVariableOp2b
/conv2d_3/bias/Regularizer/Square/ReadVariableOp/conv2d_3/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_3/kernel/Regularizer/Square/ReadVariableOp1conv2d_3/kernel/Regularizer/Square/ReadVariableOp2@
dense_2/BiasAdd/ReadVariableOpdense_2/BiasAdd/ReadVariableOp2>
dense_2/MatMul/ReadVariableOpdense_2/MatMul/ReadVariableOp2`
.dense_2/bias/Regularizer/Square/ReadVariableOp.dense_2/bias/Regularizer/Square/ReadVariableOp2d
0dense_2/kernel/Regularizer/Square/ReadVariableOp0dense_2/kernel/Regularizer/Square/ReadVariableOp2@
dense_3/BiasAdd/ReadVariableOpdense_3/BiasAdd/ReadVariableOp2>
dense_3/MatMul/ReadVariableOpdense_3/MatMul/ReadVariableOp2`
.dense_3/bias/Regularizer/Square/ReadVariableOp.dense_3/bias/Regularizer/Square/ReadVariableOp2d
0dense_3/kernel/Regularizer/Square/ReadVariableOp0dense_3/kernel/Regularizer/Square/ReadVariableOp:O K
'
_output_shapes
:?????????d
 
_user_specified_nameinputs
?
h
L__inference_up_sampling2d_1_layer_call_and_return_conditional_losses_5331516

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"      2
Constc
Const_1Const*
_output_shapes
:*
dtype0*
valueB"      2	
Const_1X
mulMulConst:output:0Const_1:output:0*
T0*
_output_shapes
:2
mul?
resize/ResizeNearestNeighborResizeNearestNeighborinputsmul:z:0*
T0*0
_output_shapes
:??????????*
half_pixel_centers(2
resize/ResizeNearestNeighbor?
IdentityIdentity-resize/ResizeNearestNeighbor:resized_images:0*
T0*0
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:??????????:X T
0
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5331346

inputs%
readvariableop_resource:@'
readvariableop_1_resource:@6
(fusedbatchnormv3_readvariableop_resource:@8
*fusedbatchnormv3_readvariableop_1_resource:@
identity??FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1t
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
:@*
dtype02
ReadVariableOpz
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes
:@*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes
:@*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes
:@*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*]
_output_shapesK
I:+???????????????????????????@:@:@:@:@:*
epsilon%o?:*
is_training( 2
FusedBatchNormV3?
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*A
_output_shapes/
-:+???????????????????????????@2

Identity?
NoOpNoOp ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*H
_input_shapes7
5:+???????????????????????????@: : : : 2B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:i e
A
_output_shapes/
-:+???????????????????????????@
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_5331175

inputs0
!batchnorm_readvariableop_resource:	?4
%batchnorm_mul_readvariableop_resource:	?2
#batchnorm_readvariableop_1_resource:	?2
#batchnorm_readvariableop_2_resource:	?
identity??batchnorm/ReadVariableOp?batchnorm/ReadVariableOp_1?batchnorm/ReadVariableOp_2?batchnorm/mul/ReadVariableOp?
batchnorm/ReadVariableOpReadVariableOp!batchnorm_readvariableop_resource*
_output_shapes	
:?*
dtype02
batchnorm/ReadVariableOpg
batchnorm/add/yConst*
_output_shapes
: *
dtype0*
valueB
 *o?:2
batchnorm/add/y?
batchnorm/addAddV2 batchnorm/ReadVariableOp:value:0batchnorm/add/y:output:0*
T0*
_output_shapes	
:?2
batchnorm/addd
batchnorm/RsqrtRsqrtbatchnorm/add:z:0*
T0*
_output_shapes	
:?2
batchnorm/Rsqrt?
batchnorm/mul/ReadVariableOpReadVariableOp%batchnorm_mul_readvariableop_resource*
_output_shapes	
:?*
dtype02
batchnorm/mul/ReadVariableOp?
batchnorm/mulMulbatchnorm/Rsqrt:y:0$batchnorm/mul/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2
batchnorm/mulw
batchnorm/mul_1Mulinputsbatchnorm/mul:z:0*
T0*(
_output_shapes
:??????????2
batchnorm/mul_1?
batchnorm/ReadVariableOp_1ReadVariableOp#batchnorm_readvariableop_1_resource*
_output_shapes	
:?*
dtype02
batchnorm/ReadVariableOp_1?
batchnorm/mul_2Mul"batchnorm/ReadVariableOp_1:value:0batchnorm/mul:z:0*
T0*
_output_shapes	
:?2
batchnorm/mul_2?
batchnorm/ReadVariableOp_2ReadVariableOp#batchnorm_readvariableop_2_resource*
_output_shapes	
:?*
dtype02
batchnorm/ReadVariableOp_2?
batchnorm/subSub"batchnorm/ReadVariableOp_2:value:0batchnorm/mul_2:z:0*
T0*
_output_shapes	
:?2
batchnorm/sub?
batchnorm/add_1AddV2batchnorm/mul_1:z:0batchnorm/sub:z:0*
T0*(
_output_shapes
:??????????2
batchnorm/add_1o
IdentityIdentitybatchnorm/add_1:z:0^NoOp*
T0*(
_output_shapes
:??????????2

Identity?
NoOpNoOp^batchnorm/ReadVariableOp^batchnorm/ReadVariableOp_1^batchnorm/ReadVariableOp_2^batchnorm/mul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*/
_input_shapes
:??????????: : : : 24
batchnorm/ReadVariableOpbatchnorm/ReadVariableOp28
batchnorm/ReadVariableOp_1batchnorm/ReadVariableOp_128
batchnorm/ReadVariableOp_2batchnorm/ReadVariableOp_22<
batchnorm/mul/ReadVariableOpbatchnorm/mul/ReadVariableOp:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
f
J__inference_up_sampling2d_layer_call_and_return_conditional_losses_5329333

inputs
identityD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2?
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2
strided_slice_
ConstConst*
_output_shapes
:*
dtype0*
valueB"      2
Const^
mulMulstrided_slice:output:0Const:output:0*
T0*
_output_shapes
:2
mul?
resize/ResizeNearestNeighborResizeNearestNeighborinputsmul:z:0*
T0*J
_output_shapes8
6:4????????????????????????????????????*
half_pixel_centers(2
resize/ResizeNearestNeighbor?
IdentityIdentity-resize/ResizeNearestNeighbor:resized_images:0*
T0*J
_output_shapes8
6:4????????????????????????????????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4????????????????????????????????????:r n
J
_output_shapes8
6:4????????????????????????????????????
 
_user_specified_nameinputs
?
d
F__inference_dropout_1_layer_call_and_return_conditional_losses_5331133

inputs

identity_1[
IdentityIdentityinputs*
T0*(
_output_shapes
:??????????2

Identityj

Identity_1IdentityIdentity:output:0*
T0*(
_output_shapes
:??????????2

Identity_1"!

identity_1Identity_1:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:??????????:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5331598

inputs&
readvariableop_resource:	?(
readvariableop_1_resource:	?7
(fusedbatchnormv3_readvariableop_resource:	?9
*fusedbatchnormv3_readvariableop_1_resource:	?
identity??AssignNewValue?AssignNewValue_1?FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1u
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp{
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes	
:?*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes	
:?*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes	
:?*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*P
_output_shapes>
<:??????????:?:?:?:?:*
epsilon%o?:*
exponential_avg_factor%
?#<2
FusedBatchNormV3?
AssignNewValueAssignVariableOp(fusedbatchnormv3_readvariableop_resourceFusedBatchNormV3:batch_mean:0 ^FusedBatchNormV3/ReadVariableOp*
_output_shapes
 *
dtype02
AssignNewValue?
AssignNewValue_1AssignVariableOp*fusedbatchnormv3_readvariableop_1_resource!FusedBatchNormV3:batch_variance:0"^FusedBatchNormV3/ReadVariableOp_1*
_output_shapes
 *
dtype02
AssignNewValue_1x
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*0
_output_shapes
:??????????2

Identity?
NoOpNoOp^AssignNewValue^AssignNewValue_1 ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*7
_input_shapes&
$:??????????: : : : 2 
AssignNewValueAssignNewValue2$
AssignNewValue_1AssignNewValue_12B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:X T
0
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
__inference_loss_fn_5_5331771G
8conv2d_3_bias_regularizer_square_readvariableop_resource:	?
identity??/conv2d_3/bias/Regularizer/Square/ReadVariableOp?
/conv2d_3/bias/Regularizer/Square/ReadVariableOpReadVariableOp8conv2d_3_bias_regularizer_square_readvariableop_resource*
_output_shapes	
:?*
dtype021
/conv2d_3/bias/Regularizer/Square/ReadVariableOp?
 conv2d_3/bias/Regularizer/SquareSquare7conv2d_3/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes	
:?2"
 conv2d_3/bias/Regularizer/Square?
conv2d_3/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_3/bias/Regularizer/Const?
conv2d_3/bias/Regularizer/SumSum$conv2d_3/bias/Regularizer/Square:y:0(conv2d_3/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/Sum?
conv2d_3/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_3/bias/Regularizer/mul/x?
conv2d_3/bias/Regularizer/mulMul(conv2d_3/bias/Regularizer/mul/x:output:0&conv2d_3/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_3/bias/Regularizer/mulk
IdentityIdentity!conv2d_3/bias/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 2

Identity?
NoOpNoOp0^conv2d_3/bias/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2b
/conv2d_3/bias/Regularizer/Square/ReadVariableOp/conv2d_3/bias/Regularizer/Square/ReadVariableOp
?
?
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5331382

inputs%
readvariableop_resource:@'
readvariableop_1_resource:@6
(fusedbatchnormv3_readvariableop_resource:@8
*fusedbatchnormv3_readvariableop_1_resource:@
identity??FusedBatchNormV3/ReadVariableOp?!FusedBatchNormV3/ReadVariableOp_1?ReadVariableOp?ReadVariableOp_1t
ReadVariableOpReadVariableOpreadvariableop_resource*
_output_shapes
:@*
dtype02
ReadVariableOpz
ReadVariableOp_1ReadVariableOpreadvariableop_1_resource*
_output_shapes
:@*
dtype02
ReadVariableOp_1?
FusedBatchNormV3/ReadVariableOpReadVariableOp(fusedbatchnormv3_readvariableop_resource*
_output_shapes
:@*
dtype02!
FusedBatchNormV3/ReadVariableOp?
!FusedBatchNormV3/ReadVariableOp_1ReadVariableOp*fusedbatchnormv3_readvariableop_1_resource*
_output_shapes
:@*
dtype02#
!FusedBatchNormV3/ReadVariableOp_1?
FusedBatchNormV3FusedBatchNormV3inputsReadVariableOp:value:0ReadVariableOp_1:value:0'FusedBatchNormV3/ReadVariableOp:value:0)FusedBatchNormV3/ReadVariableOp_1:value:0*
T0*
U0*K
_output_shapes9
7:?????????@:@:@:@:@:*
epsilon%o?:*
is_training( 2
FusedBatchNormV3w
IdentityIdentityFusedBatchNormV3:y:0^NoOp*
T0*/
_output_shapes
:?????????@2

Identity?
NoOpNoOp ^FusedBatchNormV3/ReadVariableOp"^FusedBatchNormV3/ReadVariableOp_1^ReadVariableOp^ReadVariableOp_1*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*6
_input_shapes%
#:?????????@: : : : 2B
FusedBatchNormV3/ReadVariableOpFusedBatchNormV3/ReadVariableOp2F
!FusedBatchNormV3/ReadVariableOp_1!FusedBatchNormV3/ReadVariableOp_12 
ReadVariableOpReadVariableOp2$
ReadVariableOp_1ReadVariableOp_1:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
?
?
__inference_loss_fn_2_5331738T
:conv2d_2_kernel_regularizer_square_readvariableop_resource:@@
identity??1conv2d_2/kernel/Regularizer/Square/ReadVariableOp?
1conv2d_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOp:conv2d_2_kernel_regularizer_square_readvariableop_resource*&
_output_shapes
:@@*
dtype023
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_2/kernel/Regularizer/SquareSquare9conv2d_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*&
_output_shapes
:@@2$
"conv2d_2/kernel/Regularizer/Square?
!conv2d_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_2/kernel/Regularizer/Const?
conv2d_2/kernel/Regularizer/SumSum&conv2d_2/kernel/Regularizer/Square:y:0*conv2d_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/Sum?
!conv2d_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_2/kernel/Regularizer/mul/x?
conv2d_2/kernel/Regularizer/mulMul*conv2d_2/kernel/Regularizer/mul/x:output:0(conv2d_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/mulm
IdentityIdentity#conv2d_2/kernel/Regularizer/mul:z:0^NoOp*
T0*
_output_shapes
: 2

Identity?
NoOpNoOp2^conv2d_2/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*
_input_shapes
: 2f
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp1conv2d_2/kernel/Regularizer/Square/ReadVariableOp
?
G
+__inference_reshape_1_layer_call_fn_5331254

inputs
identity?
PartitionedCallPartitionedCallinputs*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????@* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_reshape_1_layer_call_and_return_conditional_losses_53297072
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:?????????@2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*'
_input_shapes
:??????????:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
?
E__inference_conv2d_2_layer_call_and_return_conditional_losses_5331289

inputs8
conv2d_readvariableop_resource:@@-
biasadd_readvariableop_resource:@
identity??BiasAdd/ReadVariableOp?Conv2D/ReadVariableOp?/conv2d_2/bias/Regularizer/Square/ReadVariableOp?1conv2d_2/kernel/Regularizer/Square/ReadVariableOp?
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:@@*
dtype02
Conv2D/ReadVariableOp?
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????@*
paddingSAME*
strides
2
Conv2D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????@2	
BiasAdd]
EluEluBiasAdd:output:0*
T0*/
_output_shapes
:?????????@2
Elu?
1conv2d_2/kernel/Regularizer/Square/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:@@*
dtype023
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp?
"conv2d_2/kernel/Regularizer/SquareSquare9conv2d_2/kernel/Regularizer/Square/ReadVariableOp:value:0*
T0*&
_output_shapes
:@@2$
"conv2d_2/kernel/Regularizer/Square?
!conv2d_2/kernel/Regularizer/ConstConst*
_output_shapes
:*
dtype0*%
valueB"             2#
!conv2d_2/kernel/Regularizer/Const?
conv2d_2/kernel/Regularizer/SumSum&conv2d_2/kernel/Regularizer/Square:y:0*conv2d_2/kernel/Regularizer/Const:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/Sum?
!conv2d_2/kernel/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2#
!conv2d_2/kernel/Regularizer/mul/x?
conv2d_2/kernel/Regularizer/mulMul*conv2d_2/kernel/Regularizer/mul/x:output:0(conv2d_2/kernel/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2!
conv2d_2/kernel/Regularizer/mul?
/conv2d_2/bias/Regularizer/Square/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:@*
dtype021
/conv2d_2/bias/Regularizer/Square/ReadVariableOp?
 conv2d_2/bias/Regularizer/SquareSquare7conv2d_2/bias/Regularizer/Square/ReadVariableOp:value:0*
T0*
_output_shapes
:@2"
 conv2d_2/bias/Regularizer/Square?
conv2d_2/bias/Regularizer/ConstConst*
_output_shapes
:*
dtype0*
valueB: 2!
conv2d_2/bias/Regularizer/Const?
conv2d_2/bias/Regularizer/SumSum$conv2d_2/bias/Regularizer/Square:y:0(conv2d_2/bias/Regularizer/Const:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/Sum?
conv2d_2/bias/Regularizer/mul/xConst*
_output_shapes
: *
dtype0*
valueB
 *    2!
conv2d_2/bias/Regularizer/mul/x?
conv2d_2/bias/Regularizer/mulMul(conv2d_2/bias/Regularizer/mul/x:output:0&conv2d_2/bias/Regularizer/Sum:output:0*
T0*
_output_shapes
: 2
conv2d_2/bias/Regularizer/mult
IdentityIdentityElu:activations:0^NoOp*
T0*/
_output_shapes
:?????????@2

Identity?
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp0^conv2d_2/bias/Regularizer/Square/ReadVariableOp2^conv2d_2/kernel/Regularizer/Square/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????@: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp2b
/conv2d_2/bias/Regularizer/Square/ReadVariableOp/conv2d_2/bias/Regularizer/Square/ReadVariableOp2f
1conv2d_2/kernel/Regularizer/Square/ReadVariableOp1conv2d_2/kernel/Regularizer/Square/ReadVariableOp:W S
/
_output_shapes
:?????????@
 
_user_specified_nameinputs
?
f
J__inference_up_sampling2d_layer_call_and_return_conditional_losses_5331310

inputs
identityD
ShapeShapeinputs*
T0*
_output_shapes
:2
Shapet
strided_slice/stackConst*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stackx
strided_slice/stack_1Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_1x
strided_slice/stack_2Const*
_output_shapes
:*
dtype0*
valueB:2
strided_slice/stack_2?
strided_sliceStridedSliceShape:output:0strided_slice/stack:output:0strided_slice/stack_1:output:0strided_slice/stack_2:output:0*
Index0*
T0*
_output_shapes
:2
strided_slice_
ConstConst*
_output_shapes
:*
dtype0*
valueB"      2
Const^
mulMulstrided_slice:output:0Const:output:0*
T0*
_output_shapes
:2
mul?
resize/ResizeNearestNeighborResizeNearestNeighborinputsmul:z:0*
T0*J
_output_shapes8
6:4????????????????????????????????????*
half_pixel_centers(2
resize/ResizeNearestNeighbor?
IdentityIdentity-resize/ResizeNearestNeighbor:resized_images:0*
T0*J
_output_shapes8
6:4????????????????????????????????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*I
_input_shapes8
6:4????????????????????????????????????:r n
J
_output_shapes8
6:4????????????????????????????????????
 
_user_specified_nameinputs"?L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*?
serving_default?
G
dense_2_input6
serving_default_dense_2_input:0?????????d;
dense_30
StatefulPartitionedCall:0?????????2tensorflow/serving/predict:??
?
layer_with_weights-0
layer-0
layer-1
layer_with_weights-1
layer-2
layer-3
layer_with_weights-2
layer-4
layer-5
layer_with_weights-3
layer-6
layer_with_weights-4
layer-7
	layer-8

layer_with_weights-5

layer-9
layer-10
layer_with_weights-6
layer-11
trainable_variables
	variables
regularization_losses
	keras_api

signatures
+?&call_and_return_all_conditional_losses
?_default_save_signature
?__call__"
_tf_keras_sequential
?

kernel
bias
trainable_variables
	variables
regularization_losses
	keras_api
+?&call_and_return_all_conditional_losses
?__call__"
_tf_keras_layer
?
trainable_variables
	variables
regularization_losses
	keras_api
+?&call_and_return_all_conditional_losses
?__call__"
_tf_keras_layer
?
axis
	gamma
beta
moving_mean
 moving_variance
!trainable_variables
"	variables
#regularization_losses
$	keras_api
+?&call_and_return_all_conditional_losses
?__call__"
_tf_keras_layer
?
%trainable_variables
&	variables
'regularization_losses
(	keras_api
+?&call_and_return_all_conditional_losses
?__call__"
_tf_keras_layer
?

)kernel
*bias
+trainable_variables
,	variables
-regularization_losses
.	keras_api
+?&call_and_return_all_conditional_losses
?__call__"
_tf_keras_layer
?
/trainable_variables
0	variables
1regularization_losses
2	keras_api
+?&call_and_return_all_conditional_losses
?__call__"
_tf_keras_layer
?
3axis
	4gamma
5beta
6moving_mean
7moving_variance
8trainable_variables
9	variables
:regularization_losses
;	keras_api
+?&call_and_return_all_conditional_losses
?__call__"
_tf_keras_layer
?

<kernel
=bias
>trainable_variables
?	variables
@regularization_losses
A	keras_api
+?&call_and_return_all_conditional_losses
?__call__"
_tf_keras_layer
?
Btrainable_variables
C	variables
Dregularization_losses
E	keras_api
+?&call_and_return_all_conditional_losses
?__call__"
_tf_keras_layer
?
Faxis
	Ggamma
Hbeta
Imoving_mean
Jmoving_variance
Ktrainable_variables
L	variables
Mregularization_losses
N	keras_api
+?&call_and_return_all_conditional_losses
?__call__"
_tf_keras_layer
?
Otrainable_variables
P	variables
Qregularization_losses
R	keras_api
+?&call_and_return_all_conditional_losses
?__call__"
_tf_keras_layer
?

Skernel
Tbias
Utrainable_variables
V	variables
Wregularization_losses
X	keras_api
+?&call_and_return_all_conditional_losses
?__call__"
_tf_keras_layer
?
0
1
2
3
)4
*5
46
57
<8
=9
G10
H11
S12
T13"
trackable_list_wrapper
?
0
1
2
3
4
 5
)6
*7
48
59
610
711
<12
=13
G14
H15
I16
J17
S18
T19"
trackable_list_wrapper
`
?0
?1
?2
?3
?4
?5
?6
?7"
trackable_list_wrapper
?
Ymetrics
trainable_variables
	variables

Zlayers
[layer_metrics
\layer_regularization_losses
regularization_losses
]non_trainable_variables
?__call__
?_default_save_signature
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
-
?serving_default"
signature_map
!:	d?2dense_2/kernel
:?2dense_2/bias
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
0
?0
?1"
trackable_list_wrapper
?
^metrics
trainable_variables
_non_trainable_variables
	variables
`layer_metrics
alayer_regularization_losses
regularization_losses

blayers
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
cmetrics
trainable_variables
dnon_trainable_variables
	variables
elayer_metrics
flayer_regularization_losses
regularization_losses

glayers
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
*:(?2batch_normalization_3/gamma
):'?2batch_normalization_3/beta
2:0? (2!batch_normalization_3/moving_mean
6:4? (2%batch_normalization_3/moving_variance
.
0
1"
trackable_list_wrapper
<
0
1
2
 3"
trackable_list_wrapper
 "
trackable_list_wrapper
?
hmetrics
!trainable_variables
inon_trainable_variables
"	variables
jlayer_metrics
klayer_regularization_losses
#regularization_losses

llayers
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
mmetrics
%trainable_variables
nnon_trainable_variables
&	variables
olayer_metrics
player_regularization_losses
'regularization_losses

qlayers
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
):'@@2conv2d_2/kernel
:@2conv2d_2/bias
.
)0
*1"
trackable_list_wrapper
.
)0
*1"
trackable_list_wrapper
0
?0
?1"
trackable_list_wrapper
?
rmetrics
+trainable_variables
snon_trainable_variables
,	variables
tlayer_metrics
ulayer_regularization_losses
-regularization_losses

vlayers
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
wmetrics
/trainable_variables
xnon_trainable_variables
0	variables
ylayer_metrics
zlayer_regularization_losses
1regularization_losses

{layers
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
):'@2batch_normalization_4/gamma
(:&@2batch_normalization_4/beta
1:/@ (2!batch_normalization_4/moving_mean
5:3@ (2%batch_normalization_4/moving_variance
.
40
51"
trackable_list_wrapper
<
40
51
62
73"
trackable_list_wrapper
 "
trackable_list_wrapper
?
|metrics
8trainable_variables
}non_trainable_variables
9	variables
~layer_metrics
layer_regularization_losses
:regularization_losses
?layers
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
*:(@?2conv2d_3/kernel
:?2conv2d_3/bias
.
<0
=1"
trackable_list_wrapper
.
<0
=1"
trackable_list_wrapper
0
?0
?1"
trackable_list_wrapper
?
?metrics
>trainable_variables
?non_trainable_variables
?	variables
?layer_metrics
 ?layer_regularization_losses
@regularization_losses
?layers
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
Btrainable_variables
?non_trainable_variables
C	variables
?layer_metrics
 ?layer_regularization_losses
Dregularization_losses
?layers
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
*:(?2batch_normalization_5/gamma
):'?2batch_normalization_5/beta
2:0? (2!batch_normalization_5/moving_mean
6:4? (2%batch_normalization_5/moving_variance
.
G0
H1"
trackable_list_wrapper
<
G0
H1
I2
J3"
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
Ktrainable_variables
?non_trainable_variables
L	variables
?layer_metrics
 ?layer_regularization_losses
Mregularization_losses
?layers
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
?metrics
Otrainable_variables
?non_trainable_variables
P	variables
?layer_metrics
 ?layer_regularization_losses
Qregularization_losses
?layers
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
!:	?@22dense_3/kernel
:22dense_3/bias
.
S0
T1"
trackable_list_wrapper
.
S0
T1"
trackable_list_wrapper
0
?0
?1"
trackable_list_wrapper
?
?metrics
Utrainable_variables
?non_trainable_variables
V	variables
?layer_metrics
 ?layer_regularization_losses
Wregularization_losses
?layers
?__call__
+?&call_and_return_all_conditional_losses
'?"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
v
0
1
2
3
4
5
6
7
	8

9
10
11"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
J
0
 1
62
73
I4
J5"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
0
?0
?1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
0
 1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
0
?0
?1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
60
71"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
0
?0
?1"
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
.
I0
J1"
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_dict_wrapper
0
?0
?1"
trackable_list_wrapper
 "
trackable_list_wrapper
?2?
I__inference_sequential_1_layer_call_and_return_conditional_losses_5330828
I__inference_sequential_1_layer_call_and_return_conditional_losses_5330994
I__inference_sequential_1_layer_call_and_return_conditional_losses_5330484
I__inference_sequential_1_layer_call_and_return_conditional_losses_5330588?
???
FullArgSpec1
args)?&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults?
p 

 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?B?
"__inference__wrapped_model_5329155dense_2_input"?
???
FullArgSpec
args? 
varargsjargs
varkwjkwargs
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
.__inference_sequential_1_layer_call_fn_5329968
.__inference_sequential_1_layer_call_fn_5331039
.__inference_sequential_1_layer_call_fn_5331084
.__inference_sequential_1_layer_call_fn_5330380?
???
FullArgSpec1
args)?&
jself
jinputs

jtraining
jmask
varargs
 
varkw
 
defaults?
p 

 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
D__inference_dense_2_layer_call_and_return_conditional_losses_5331119?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
)__inference_dense_2_layer_call_fn_5331128?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
F__inference_dropout_1_layer_call_and_return_conditional_losses_5331133
F__inference_dropout_1_layer_call_and_return_conditional_losses_5331145?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
+__inference_dropout_1_layer_call_fn_5331150
+__inference_dropout_1_layer_call_fn_5331155?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_5331175
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_5331209?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
7__inference_batch_normalization_3_layer_call_fn_5331222
7__inference_batch_normalization_3_layer_call_fn_5331235?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
F__inference_reshape_1_layer_call_and_return_conditional_losses_5331249?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
+__inference_reshape_1_layer_call_fn_5331254?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
E__inference_conv2d_2_layer_call_and_return_conditional_losses_5331289?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
*__inference_conv2d_2_layer_call_fn_5331298?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
J__inference_up_sampling2d_layer_call_and_return_conditional_losses_5331310
J__inference_up_sampling2d_layer_call_and_return_conditional_losses_5331318?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
/__inference_up_sampling2d_layer_call_fn_5331323
/__inference_up_sampling2d_layer_call_fn_5331328?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5331346
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5331364
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5331382
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5331400?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
7__inference_batch_normalization_4_layer_call_fn_5331413
7__inference_batch_normalization_4_layer_call_fn_5331426
7__inference_batch_normalization_4_layer_call_fn_5331439
7__inference_batch_normalization_4_layer_call_fn_5331452?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
E__inference_conv2d_3_layer_call_and_return_conditional_losses_5331487?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
*__inference_conv2d_3_layer_call_fn_5331496?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
L__inference_up_sampling2d_1_layer_call_and_return_conditional_losses_5331508
L__inference_up_sampling2d_1_layer_call_and_return_conditional_losses_5331516?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
1__inference_up_sampling2d_1_layer_call_fn_5331521
1__inference_up_sampling2d_1_layer_call_fn_5331526?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5331544
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5331562
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5331580
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5331598?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
7__inference_batch_normalization_5_layer_call_fn_5331611
7__inference_batch_normalization_5_layer_call_fn_5331624
7__inference_batch_normalization_5_layer_call_fn_5331637
7__inference_batch_normalization_5_layer_call_fn_5331650?
???
FullArgSpec)
args!?
jself
jinputs

jtraining
varargs
 
varkw
 
defaults?
p 

kwonlyargs? 
kwonlydefaults? 
annotations? *
 
?2?
F__inference_flatten_1_layer_call_and_return_conditional_losses_5331656?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
+__inference_flatten_1_layer_call_fn_5331661?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
D__inference_dense_3_layer_call_and_return_conditional_losses_5331696?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
)__inference_dense_3_layer_call_fn_5331705?
???
FullArgSpec
args?
jself
jinputs
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 
?2?
__inference_loss_fn_0_5331716?
???
FullArgSpec
args? 
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *? 
?2?
__inference_loss_fn_1_5331727?
???
FullArgSpec
args? 
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *? 
?2?
__inference_loss_fn_2_5331738?
???
FullArgSpec
args? 
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *? 
?2?
__inference_loss_fn_3_5331749?
???
FullArgSpec
args? 
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *? 
?2?
__inference_loss_fn_4_5331760?
???
FullArgSpec
args? 
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *? 
?2?
__inference_loss_fn_5_5331771?
???
FullArgSpec
args? 
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *? 
?2?
__inference_loss_fn_6_5331782?
???
FullArgSpec
args? 
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *? 
?2?
__inference_loss_fn_7_5331793?
???
FullArgSpec
args? 
varargs
 
varkw
 
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *? 
?B?
%__inference_signature_wrapper_5330683dense_2_input"?
???
FullArgSpec
args? 
varargs
 
varkwjkwargs
defaults
 

kwonlyargs? 
kwonlydefaults
 
annotations? *
 ?
"__inference__wrapped_model_5329155? )*4567<=GHIJST6?3
,?)
'?$
dense_2_input?????????d
? "1?.
,
dense_3!?
dense_3?????????2?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_5331175d 4?1
*?'
!?
inputs??????????
p 
? "&?#
?
0??????????
? ?
R__inference_batch_normalization_3_layer_call_and_return_conditional_losses_5331209d 4?1
*?'
!?
inputs??????????
p
? "&?#
?
0??????????
? ?
7__inference_batch_normalization_3_layer_call_fn_5331222W 4?1
*?'
!?
inputs??????????
p 
? "????????????
7__inference_batch_normalization_3_layer_call_fn_5331235W 4?1
*?'
!?
inputs??????????
p
? "????????????
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5331346?4567M?J
C?@
:?7
inputs+???????????????????????????@
p 
? "??<
5?2
0+???????????????????????????@
? ?
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5331364?4567M?J
C?@
:?7
inputs+???????????????????????????@
p
? "??<
5?2
0+???????????????????????????@
? ?
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5331382r4567;?8
1?.
(?%
inputs?????????@
p 
? "-?*
#? 
0?????????@
? ?
R__inference_batch_normalization_4_layer_call_and_return_conditional_losses_5331400r4567;?8
1?.
(?%
inputs?????????@
p
? "-?*
#? 
0?????????@
? ?
7__inference_batch_normalization_4_layer_call_fn_5331413?4567M?J
C?@
:?7
inputs+???????????????????????????@
p 
? "2?/+???????????????????????????@?
7__inference_batch_normalization_4_layer_call_fn_5331426?4567M?J
C?@
:?7
inputs+???????????????????????????@
p
? "2?/+???????????????????????????@?
7__inference_batch_normalization_4_layer_call_fn_5331439e4567;?8
1?.
(?%
inputs?????????@
p 
? " ??????????@?
7__inference_batch_normalization_4_layer_call_fn_5331452e4567;?8
1?.
(?%
inputs?????????@
p
? " ??????????@?
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5331544?GHIJN?K
D?A
;?8
inputs,????????????????????????????
p 
? "@?=
6?3
0,????????????????????????????
? ?
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5331562?GHIJN?K
D?A
;?8
inputs,????????????????????????????
p
? "@?=
6?3
0,????????????????????????????
? ?
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5331580tGHIJ<?9
2?/
)?&
inputs??????????
p 
? ".?+
$?!
0??????????
? ?
R__inference_batch_normalization_5_layer_call_and_return_conditional_losses_5331598tGHIJ<?9
2?/
)?&
inputs??????????
p
? ".?+
$?!
0??????????
? ?
7__inference_batch_normalization_5_layer_call_fn_5331611?GHIJN?K
D?A
;?8
inputs,????????????????????????????
p 
? "3?0,?????????????????????????????
7__inference_batch_normalization_5_layer_call_fn_5331624?GHIJN?K
D?A
;?8
inputs,????????????????????????????
p
? "3?0,?????????????????????????????
7__inference_batch_normalization_5_layer_call_fn_5331637gGHIJ<?9
2?/
)?&
inputs??????????
p 
? "!????????????
7__inference_batch_normalization_5_layer_call_fn_5331650gGHIJ<?9
2?/
)?&
inputs??????????
p
? "!????????????
E__inference_conv2d_2_layer_call_and_return_conditional_losses_5331289l)*7?4
-?*
(?%
inputs?????????@
? "-?*
#? 
0?????????@
? ?
*__inference_conv2d_2_layer_call_fn_5331298_)*7?4
-?*
(?%
inputs?????????@
? " ??????????@?
E__inference_conv2d_3_layer_call_and_return_conditional_losses_5331487m<=7?4
-?*
(?%
inputs?????????@
? ".?+
$?!
0??????????
? ?
*__inference_conv2d_3_layer_call_fn_5331496`<=7?4
-?*
(?%
inputs?????????@
? "!????????????
D__inference_dense_2_layer_call_and_return_conditional_losses_5331119]/?,
%?"
 ?
inputs?????????d
? "&?#
?
0??????????
? }
)__inference_dense_2_layer_call_fn_5331128P/?,
%?"
 ?
inputs?????????d
? "????????????
D__inference_dense_3_layer_call_and_return_conditional_losses_5331696]ST0?-
&?#
!?
inputs??????????@
? "%?"
?
0?????????2
? }
)__inference_dense_3_layer_call_fn_5331705PST0?-
&?#
!?
inputs??????????@
? "??????????2?
F__inference_dropout_1_layer_call_and_return_conditional_losses_5331133^4?1
*?'
!?
inputs??????????
p 
? "&?#
?
0??????????
? ?
F__inference_dropout_1_layer_call_and_return_conditional_losses_5331145^4?1
*?'
!?
inputs??????????
p
? "&?#
?
0??????????
? ?
+__inference_dropout_1_layer_call_fn_5331150Q4?1
*?'
!?
inputs??????????
p 
? "????????????
+__inference_dropout_1_layer_call_fn_5331155Q4?1
*?'
!?
inputs??????????
p
? "????????????
F__inference_flatten_1_layer_call_and_return_conditional_losses_5331656b8?5
.?+
)?&
inputs??????????
? "&?#
?
0??????????@
? ?
+__inference_flatten_1_layer_call_fn_5331661U8?5
.?+
)?&
inputs??????????
? "???????????@<
__inference_loss_fn_0_5331716?

? 
? "? <
__inference_loss_fn_1_5331727?

? 
? "? <
__inference_loss_fn_2_5331738)?

? 
? "? <
__inference_loss_fn_3_5331749*?

? 
? "? <
__inference_loss_fn_4_5331760<?

? 
? "? <
__inference_loss_fn_5_5331771=?

? 
? "? <
__inference_loss_fn_6_5331782S?

? 
? "? <
__inference_loss_fn_7_5331793T?

? 
? "? ?
F__inference_reshape_1_layer_call_and_return_conditional_losses_5331249a0?-
&?#
!?
inputs??????????
? "-?*
#? 
0?????????@
? ?
+__inference_reshape_1_layer_call_fn_5331254T0?-
&?#
!?
inputs??????????
? " ??????????@?
I__inference_sequential_1_layer_call_and_return_conditional_losses_5330484} )*4567<=GHIJST>?;
4?1
'?$
dense_2_input?????????d
p 

 
? "%?"
?
0?????????2
? ?
I__inference_sequential_1_layer_call_and_return_conditional_losses_5330588} )*4567<=GHIJST>?;
4?1
'?$
dense_2_input?????????d
p

 
? "%?"
?
0?????????2
? ?
I__inference_sequential_1_layer_call_and_return_conditional_losses_5330828v )*4567<=GHIJST7?4
-?*
 ?
inputs?????????d
p 

 
? "%?"
?
0?????????2
? ?
I__inference_sequential_1_layer_call_and_return_conditional_losses_5330994v )*4567<=GHIJST7?4
-?*
 ?
inputs?????????d
p

 
? "%?"
?
0?????????2
? ?
.__inference_sequential_1_layer_call_fn_5329968p )*4567<=GHIJST>?;
4?1
'?$
dense_2_input?????????d
p 

 
? "??????????2?
.__inference_sequential_1_layer_call_fn_5330380p )*4567<=GHIJST>?;
4?1
'?$
dense_2_input?????????d
p

 
? "??????????2?
.__inference_sequential_1_layer_call_fn_5331039i )*4567<=GHIJST7?4
-?*
 ?
inputs?????????d
p 

 
? "??????????2?
.__inference_sequential_1_layer_call_fn_5331084i )*4567<=GHIJST7?4
-?*
 ?
inputs?????????d
p

 
? "??????????2?
%__inference_signature_wrapper_5330683? )*4567<=GHIJSTG?D
? 
=?:
8
dense_2_input'?$
dense_2_input?????????d"1?.
,
dense_3!?
dense_3?????????2?
L__inference_up_sampling2d_1_layer_call_and_return_conditional_losses_5331508?R?O
H?E
C?@
inputs4????????????????????????????????????
? "H?E
>?;
04????????????????????????????????????
? ?
L__inference_up_sampling2d_1_layer_call_and_return_conditional_losses_5331516j8?5
.?+
)?&
inputs??????????
? ".?+
$?!
0??????????
? ?
1__inference_up_sampling2d_1_layer_call_fn_5331521?R?O
H?E
C?@
inputs4????????????????????????????????????
? ";?84?????????????????????????????????????
1__inference_up_sampling2d_1_layer_call_fn_5331526]8?5
.?+
)?&
inputs??????????
? "!????????????
J__inference_up_sampling2d_layer_call_and_return_conditional_losses_5331310?R?O
H?E
C?@
inputs4????????????????????????????????????
? "H?E
>?;
04????????????????????????????????????
? ?
J__inference_up_sampling2d_layer_call_and_return_conditional_losses_5331318h7?4
-?*
(?%
inputs?????????@
? "-?*
#? 
0?????????@
? ?
/__inference_up_sampling2d_layer_call_fn_5331323?R?O
H?E
C?@
inputs4????????????????????????????????????
? ";?84?????????????????????????????????????
/__inference_up_sampling2d_layer_call_fn_5331328[7?4
-?*
(?%
inputs?????????@
? " ??????????@