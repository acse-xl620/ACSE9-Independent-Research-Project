܁
??
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
?
MaxPool

input"T
output"T"
Ttype0:
2	"
ksize	list(int)(0"
strides	list(int)(0",
paddingstring:
SAMEVALIDEXPLICIT""
explicit_paddings	list(int)
 ":
data_formatstringNHWC:
NHWCNCHWNCHW_VECT_C
e
MergeV2Checkpoints
checkpoint_prefixes
destination_prefix"
delete_old_dirsbool(?
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
o
	RestoreV2

prefix
tensor_names
shape_and_slices
tensors2dtypes"
dtypes
list(type)(0?
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
H
ShardedFilename
basename	
shard

num_shards
filename
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
N

StringJoin
inputs*N

output"
Nint(0"
	separatorstring 
?
VarHandleOp
resource"
	containerstring "
shared_namestring "
dtypetype"
shapeshape"#
allowed_deviceslist(string)
 ?"serve*2.6.02v2.6.0-0-g919f693420e8??
?
conv2d_14/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:*!
shared_nameconv2d_14/kernel
}
$conv2d_14/kernel/Read/ReadVariableOpReadVariableOpconv2d_14/kernel*&
_output_shapes
:*
dtype0
t
conv2d_14/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameconv2d_14/bias
m
"conv2d_14/bias/Read/ReadVariableOpReadVariableOpconv2d_14/bias*
_output_shapes
:*
dtype0
?
conv2d_15/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:*!
shared_nameconv2d_15/kernel
}
$conv2d_15/kernel/Read/ReadVariableOpReadVariableOpconv2d_15/kernel*&
_output_shapes
:*
dtype0
t
conv2d_15/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameconv2d_15/bias
m
"conv2d_15/bias/Read/ReadVariableOpReadVariableOpconv2d_15/bias*
_output_shapes
:*
dtype0
?
conv2d_16/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:*!
shared_nameconv2d_16/kernel
}
$conv2d_16/kernel/Read/ReadVariableOpReadVariableOpconv2d_16/kernel*&
_output_shapes
:*
dtype0
t
conv2d_16/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_nameconv2d_16/bias
m
"conv2d_16/bias/Read/ReadVariableOpReadVariableOpconv2d_16/bias*
_output_shapes
:*
dtype0
y
dense_4/kernelVarHandleOp*
_output_shapes
: *
dtype0*
shape:	?*
shared_namedense_4/kernel
r
"dense_4/kernel/Read/ReadVariableOpReadVariableOpdense_4/kernel*
_output_shapes
:	?*
dtype0
p
dense_4/biasVarHandleOp*
_output_shapes
: *
dtype0*
shape:*
shared_namedense_4/bias
i
 dense_4/bias/Read/ReadVariableOpReadVariableOpdense_4/bias*
_output_shapes
:*
dtype0

NoOpNoOp
?
ConstConst"/device:CPU:0*
_output_shapes
: *
dtype0*?
value?B? B?
?
layer_with_weights-0
layer-0
layer-1
layer_with_weights-1
layer-2
layer-3
layer_with_weights-2
layer-4
layer-5
layer-6
layer_with_weights-3
layer-7
	regularization_losses

trainable_variables
	variables
	keras_api

signatures
h

kernel
bias
regularization_losses
trainable_variables
	variables
	keras_api
R
regularization_losses
trainable_variables
	variables
	keras_api
h

kernel
bias
regularization_losses
trainable_variables
	variables
	keras_api
R
regularization_losses
trainable_variables
 	variables
!	keras_api
h

"kernel
#bias
$regularization_losses
%trainable_variables
&	variables
'	keras_api
R
(regularization_losses
)trainable_variables
*	variables
+	keras_api
R
,regularization_losses
-trainable_variables
.	variables
/	keras_api
h

0kernel
1bias
2regularization_losses
3trainable_variables
4	variables
5	keras_api
 
8
0
1
2
3
"4
#5
06
17
8
0
1
2
3
"4
#5
06
17
?
	regularization_losses
6layer_regularization_losses

trainable_variables
	variables

7layers
8non_trainable_variables
9layer_metrics
:metrics
 
\Z
VARIABLE_VALUEconv2d_14/kernel6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUE
XV
VARIABLE_VALUEconv2d_14/bias4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUE
 

0
1

0
1
?
regularization_losses
;layer_regularization_losses
trainable_variables
	variables

<layers
=non_trainable_variables
>layer_metrics
?metrics
 
 
 
?
regularization_losses
@layer_regularization_losses
trainable_variables
	variables

Alayers
Bnon_trainable_variables
Clayer_metrics
Dmetrics
\Z
VARIABLE_VALUEconv2d_15/kernel6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUE
XV
VARIABLE_VALUEconv2d_15/bias4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUE
 

0
1

0
1
?
regularization_losses
Elayer_regularization_losses
trainable_variables
	variables

Flayers
Gnon_trainable_variables
Hlayer_metrics
Imetrics
 
 
 
?
regularization_losses
Jlayer_regularization_losses
trainable_variables
 	variables

Klayers
Lnon_trainable_variables
Mlayer_metrics
Nmetrics
\Z
VARIABLE_VALUEconv2d_16/kernel6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUE
XV
VARIABLE_VALUEconv2d_16/bias4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUE
 

"0
#1

"0
#1
?
$regularization_losses
Olayer_regularization_losses
%trainable_variables
&	variables

Players
Qnon_trainable_variables
Rlayer_metrics
Smetrics
 
 
 
?
(regularization_losses
Tlayer_regularization_losses
)trainable_variables
*	variables

Ulayers
Vnon_trainable_variables
Wlayer_metrics
Xmetrics
 
 
 
?
,regularization_losses
Ylayer_regularization_losses
-trainable_variables
.	variables

Zlayers
[non_trainable_variables
\layer_metrics
]metrics
ZX
VARIABLE_VALUEdense_4/kernel6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUE
VT
VARIABLE_VALUEdense_4/bias4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUE
 

00
11

00
11
?
2regularization_losses
^layer_regularization_losses
3trainable_variables
4	variables

_layers
`non_trainable_variables
alayer_metrics
bmetrics
 
8
0
1
2
3
4
5
6
7
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
 
 
 
 
 
?
serving_default_conv2d_14_inputPlaceholder*/
_output_shapes
:?????????00*
dtype0*$
shape:?????????00
?
StatefulPartitionedCallStatefulPartitionedCallserving_default_conv2d_14_inputconv2d_14/kernelconv2d_14/biasconv2d_15/kernelconv2d_15/biasconv2d_16/kernelconv2d_16/biasdense_4/kerneldense_4/bias*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????**
_read_only_resource_inputs

*0
config_proto 

CPU

GPU2*0J 8? *.
f)R'
%__inference_signature_wrapper_8269913
O
saver_filenamePlaceholder*
_output_shapes
: *
dtype0*
shape: 
?
StatefulPartitionedCall_1StatefulPartitionedCallsaver_filename$conv2d_14/kernel/Read/ReadVariableOp"conv2d_14/bias/Read/ReadVariableOp$conv2d_15/kernel/Read/ReadVariableOp"conv2d_15/bias/Read/ReadVariableOp$conv2d_16/kernel/Read/ReadVariableOp"conv2d_16/bias/Read/ReadVariableOp"dense_4/kernel/Read/ReadVariableOp dense_4/bias/Read/ReadVariableOpConst*
Tin
2
*
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
 __inference__traced_save_8270224
?
StatefulPartitionedCall_2StatefulPartitionedCallsaver_filenameconv2d_14/kernelconv2d_14/biasconv2d_15/kernelconv2d_15/biasconv2d_16/kernelconv2d_16/biasdense_4/kerneldense_4/bias*
Tin
2	*
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
#__inference__traced_restore_8270258??
?

?
.__inference_sequential_4_layer_call_fn_8270006

inputs!
unknown:
	unknown_0:#
	unknown_1:
	unknown_2:#
	unknown_3:
	unknown_4:
	unknown_5:	?
	unknown_6:
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????**
_read_only_resource_inputs

*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_82696632
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:?????????00: : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:?????????00
 
_user_specified_nameinputs
?
h
L__inference_max_pooling2d_6_layer_call_and_return_conditional_losses_8270057

inputs
identity?
MaxPoolMaxPoolinputs*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2	
MaxPooll
IdentityIdentityMaxPool:output:0*
T0*/
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????00:W S
/
_output_shapes
:?????????00
 
_user_specified_nameinputs
?=
?
"__inference__wrapped_model_8269496
conv2d_14_inputO
5sequential_4_conv2d_14_conv2d_readvariableop_resource:D
6sequential_4_conv2d_14_biasadd_readvariableop_resource:O
5sequential_4_conv2d_15_conv2d_readvariableop_resource:D
6sequential_4_conv2d_15_biasadd_readvariableop_resource:O
5sequential_4_conv2d_16_conv2d_readvariableop_resource:D
6sequential_4_conv2d_16_biasadd_readvariableop_resource:F
3sequential_4_dense_4_matmul_readvariableop_resource:	?B
4sequential_4_dense_4_biasadd_readvariableop_resource:
identity??-sequential_4/conv2d_14/BiasAdd/ReadVariableOp?,sequential_4/conv2d_14/Conv2D/ReadVariableOp?-sequential_4/conv2d_15/BiasAdd/ReadVariableOp?,sequential_4/conv2d_15/Conv2D/ReadVariableOp?-sequential_4/conv2d_16/BiasAdd/ReadVariableOp?,sequential_4/conv2d_16/Conv2D/ReadVariableOp?+sequential_4/dense_4/BiasAdd/ReadVariableOp?*sequential_4/dense_4/MatMul/ReadVariableOp?
,sequential_4/conv2d_14/Conv2D/ReadVariableOpReadVariableOp5sequential_4_conv2d_14_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02.
,sequential_4/conv2d_14/Conv2D/ReadVariableOp?
sequential_4/conv2d_14/Conv2DConv2Dconv2d_14_input4sequential_4/conv2d_14/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????00*
paddingSAME*
strides
2
sequential_4/conv2d_14/Conv2D?
-sequential_4/conv2d_14/BiasAdd/ReadVariableOpReadVariableOp6sequential_4_conv2d_14_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02/
-sequential_4/conv2d_14/BiasAdd/ReadVariableOp?
sequential_4/conv2d_14/BiasAddBiasAdd&sequential_4/conv2d_14/Conv2D:output:05sequential_4/conv2d_14/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????002 
sequential_4/conv2d_14/BiasAdd?
sequential_4/conv2d_14/EluElu'sequential_4/conv2d_14/BiasAdd:output:0*
T0*/
_output_shapes
:?????????002
sequential_4/conv2d_14/Elu?
$sequential_4/max_pooling2d_6/MaxPoolMaxPool(sequential_4/conv2d_14/Elu:activations:0*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2&
$sequential_4/max_pooling2d_6/MaxPool?
,sequential_4/conv2d_15/Conv2D/ReadVariableOpReadVariableOp5sequential_4_conv2d_15_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02.
,sequential_4/conv2d_15/Conv2D/ReadVariableOp?
sequential_4/conv2d_15/Conv2DConv2D-sequential_4/max_pooling2d_6/MaxPool:output:04sequential_4/conv2d_15/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????*
paddingSAME*
strides
2
sequential_4/conv2d_15/Conv2D?
-sequential_4/conv2d_15/BiasAdd/ReadVariableOpReadVariableOp6sequential_4_conv2d_15_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02/
-sequential_4/conv2d_15/BiasAdd/ReadVariableOp?
sequential_4/conv2d_15/BiasAddBiasAdd&sequential_4/conv2d_15/Conv2D:output:05sequential_4/conv2d_15/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????2 
sequential_4/conv2d_15/BiasAdd?
sequential_4/conv2d_15/EluElu'sequential_4/conv2d_15/BiasAdd:output:0*
T0*/
_output_shapes
:?????????2
sequential_4/conv2d_15/Elu?
$sequential_4/max_pooling2d_7/MaxPoolMaxPool(sequential_4/conv2d_15/Elu:activations:0*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2&
$sequential_4/max_pooling2d_7/MaxPool?
,sequential_4/conv2d_16/Conv2D/ReadVariableOpReadVariableOp5sequential_4_conv2d_16_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02.
,sequential_4/conv2d_16/Conv2D/ReadVariableOp?
sequential_4/conv2d_16/Conv2DConv2D-sequential_4/max_pooling2d_7/MaxPool:output:04sequential_4/conv2d_16/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????*
paddingSAME*
strides
2
sequential_4/conv2d_16/Conv2D?
-sequential_4/conv2d_16/BiasAdd/ReadVariableOpReadVariableOp6sequential_4_conv2d_16_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02/
-sequential_4/conv2d_16/BiasAdd/ReadVariableOp?
sequential_4/conv2d_16/BiasAddBiasAdd&sequential_4/conv2d_16/Conv2D:output:05sequential_4/conv2d_16/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????2 
sequential_4/conv2d_16/BiasAdd?
sequential_4/conv2d_16/EluElu'sequential_4/conv2d_16/BiasAdd:output:0*
T0*/
_output_shapes
:?????????2
sequential_4/conv2d_16/Elu?
$sequential_4/max_pooling2d_8/MaxPoolMaxPool(sequential_4/conv2d_16/Elu:activations:0*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2&
$sequential_4/max_pooling2d_8/MaxPool?
sequential_4/flatten_2/ConstConst*
_output_shapes
:*
dtype0*
valueB"????   2
sequential_4/flatten_2/Const?
sequential_4/flatten_2/ReshapeReshape-sequential_4/max_pooling2d_8/MaxPool:output:0%sequential_4/flatten_2/Const:output:0*
T0*(
_output_shapes
:??????????2 
sequential_4/flatten_2/Reshape?
*sequential_4/dense_4/MatMul/ReadVariableOpReadVariableOp3sequential_4_dense_4_matmul_readvariableop_resource*
_output_shapes
:	?*
dtype02,
*sequential_4/dense_4/MatMul/ReadVariableOp?
sequential_4/dense_4/MatMulMatMul'sequential_4/flatten_2/Reshape:output:02sequential_4/dense_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
sequential_4/dense_4/MatMul?
+sequential_4/dense_4/BiasAdd/ReadVariableOpReadVariableOp4sequential_4_dense_4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02-
+sequential_4/dense_4/BiasAdd/ReadVariableOp?
sequential_4/dense_4/BiasAddBiasAdd%sequential_4/dense_4/MatMul:product:03sequential_4/dense_4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
sequential_4/dense_4/BiasAdd?
IdentityIdentity%sequential_4/dense_4/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identity?
NoOpNoOp.^sequential_4/conv2d_14/BiasAdd/ReadVariableOp-^sequential_4/conv2d_14/Conv2D/ReadVariableOp.^sequential_4/conv2d_15/BiasAdd/ReadVariableOp-^sequential_4/conv2d_15/Conv2D/ReadVariableOp.^sequential_4/conv2d_16/BiasAdd/ReadVariableOp-^sequential_4/conv2d_16/Conv2D/ReadVariableOp,^sequential_4/dense_4/BiasAdd/ReadVariableOp+^sequential_4/dense_4/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:?????????00: : : : : : : : 2^
-sequential_4/conv2d_14/BiasAdd/ReadVariableOp-sequential_4/conv2d_14/BiasAdd/ReadVariableOp2\
,sequential_4/conv2d_14/Conv2D/ReadVariableOp,sequential_4/conv2d_14/Conv2D/ReadVariableOp2^
-sequential_4/conv2d_15/BiasAdd/ReadVariableOp-sequential_4/conv2d_15/BiasAdd/ReadVariableOp2\
,sequential_4/conv2d_15/Conv2D/ReadVariableOp,sequential_4/conv2d_15/Conv2D/ReadVariableOp2^
-sequential_4/conv2d_16/BiasAdd/ReadVariableOp-sequential_4/conv2d_16/BiasAdd/ReadVariableOp2\
,sequential_4/conv2d_16/Conv2D/ReadVariableOp,sequential_4/conv2d_16/Conv2D/ReadVariableOp2Z
+sequential_4/dense_4/BiasAdd/ReadVariableOp+sequential_4/dense_4/BiasAdd/ReadVariableOp2X
*sequential_4/dense_4/MatMul/ReadVariableOp*sequential_4/dense_4/MatMul/ReadVariableOp:` \
/
_output_shapes
:?????????00
)
_user_specified_nameconv2d_14_input
?
b
F__inference_flatten_2_layer_call_and_return_conditional_losses_8269644

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"????   2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:??????????2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?
?
 __inference__traced_save_8270224
file_prefix/
+savev2_conv2d_14_kernel_read_readvariableop-
)savev2_conv2d_14_bias_read_readvariableop/
+savev2_conv2d_15_kernel_read_readvariableop-
)savev2_conv2d_15_bias_read_readvariableop/
+savev2_conv2d_16_kernel_read_readvariableop-
)savev2_conv2d_16_bias_read_readvariableop-
)savev2_dense_4_kernel_read_readvariableop+
'savev2_dense_4_bias_read_readvariableop
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
ShardedFilename?
SaveV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:	*
dtype0*?
value?B?	B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
SaveV2/tensor_names?
SaveV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:	*
dtype0*%
valueB	B B B B B B B B B 2
SaveV2/shape_and_slices?
SaveV2SaveV2ShardedFilename:filename:0SaveV2/tensor_names:output:0 SaveV2/shape_and_slices:output:0+savev2_conv2d_14_kernel_read_readvariableop)savev2_conv2d_14_bias_read_readvariableop+savev2_conv2d_15_kernel_read_readvariableop)savev2_conv2d_15_bias_read_readvariableop+savev2_conv2d_16_kernel_read_readvariableop)savev2_conv2d_16_bias_read_readvariableop)savev2_dense_4_kernel_read_readvariableop'savev2_dense_4_bias_read_readvariableopsavev2_const"/device:CPU:0*
_output_shapes
 *
dtypes
2	2
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

identity_1Identity_1:output:0*p
_input_shapes_
]: :::::::	?:: 2(
MergeV2CheckpointsMergeV2Checkpoints:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix:,(
&
_output_shapes
:: 

_output_shapes
::,(
&
_output_shapes
:: 

_output_shapes
::,(
&
_output_shapes
:: 

_output_shapes
::%!

_output_shapes
:	?: 

_output_shapes
::	

_output_shapes
: 
?

?
.__inference_sequential_4_layer_call_fn_8269682
conv2d_14_input!
unknown:
	unknown_0:#
	unknown_1:
	unknown_2:#
	unknown_3:
	unknown_4:
	unknown_5:	?
	unknown_6:
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallconv2d_14_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????**
_read_only_resource_inputs

*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_82696632
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:?????????00: : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:` \
/
_output_shapes
:?????????00
)
_user_specified_nameconv2d_14_input
?
?
F__inference_conv2d_14_layer_call_and_return_conditional_losses_8270038

inputs8
conv2d_readvariableop_resource:-
biasadd_readvariableop_resource:
identity??BiasAdd/ReadVariableOp?Conv2D/ReadVariableOp?
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
Conv2D/ReadVariableOp?
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????00*
paddingSAME*
strides
2
Conv2D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????002	
BiasAdd]
EluEluBiasAdd:output:0*
T0*/
_output_shapes
:?????????002
Elut
IdentityIdentityElu:activations:0^NoOp*
T0*/
_output_shapes
:?????????002

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????00: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:?????????00
 
_user_specified_nameinputs
?
h
L__inference_max_pooling2d_7_layer_call_and_return_conditional_losses_8269613

inputs
identity?
MaxPoolMaxPoolinputs*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2	
MaxPooll
IdentityIdentityMaxPool:output:0*
T0*/
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?
h
L__inference_max_pooling2d_7_layer_call_and_return_conditional_losses_8269527

inputs
identity?
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4????????????????????????????????????*
ksize
*
paddingSAME*
strides
2	
MaxPool?
IdentityIdentityMaxPool:output:0*
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
F__inference_conv2d_16_layer_call_and_return_conditional_losses_8270118

inputs8
conv2d_readvariableop_resource:-
biasadd_readvariableop_resource:
identity??BiasAdd/ReadVariableOp?Conv2D/ReadVariableOp?
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
Conv2D/ReadVariableOp?
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????*
paddingSAME*
strides
2
Conv2D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????2	
BiasAdd]
EluEluBiasAdd:output:0*
T0*/
_output_shapes
:?????????2
Elut
IdentityIdentityElu:activations:0^NoOp*
T0*/
_output_shapes
:?????????2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?%
?
I__inference_sequential_4_layer_call_and_return_conditional_losses_8269663

inputs+
conv2d_14_8269581:
conv2d_14_8269583:+
conv2d_15_8269604:
conv2d_15_8269606:+
conv2d_16_8269627:
conv2d_16_8269629:"
dense_4_8269657:	?
dense_4_8269659:
identity??!conv2d_14/StatefulPartitionedCall?!conv2d_15/StatefulPartitionedCall?!conv2d_16/StatefulPartitionedCall?dense_4/StatefulPartitionedCall?
!conv2d_14/StatefulPartitionedCallStatefulPartitionedCallinputsconv2d_14_8269581conv2d_14_8269583*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????00*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_14_layer_call_and_return_conditional_losses_82695802#
!conv2d_14/StatefulPartitionedCall?
max_pooling2d_6/PartitionedCallPartitionedCall*conv2d_14/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_6_layer_call_and_return_conditional_losses_82695902!
max_pooling2d_6/PartitionedCall?
!conv2d_15/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_6/PartitionedCall:output:0conv2d_15_8269604conv2d_15_8269606*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_15_layer_call_and_return_conditional_losses_82696032#
!conv2d_15/StatefulPartitionedCall?
max_pooling2d_7/PartitionedCallPartitionedCall*conv2d_15/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_7_layer_call_and_return_conditional_losses_82696132!
max_pooling2d_7/PartitionedCall?
!conv2d_16/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_7/PartitionedCall:output:0conv2d_16_8269627conv2d_16_8269629*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_16_layer_call_and_return_conditional_losses_82696262#
!conv2d_16/StatefulPartitionedCall?
max_pooling2d_8/PartitionedCallPartitionedCall*conv2d_16/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_8_layer_call_and_return_conditional_losses_82696362!
max_pooling2d_8/PartitionedCall?
flatten_2/PartitionedCallPartitionedCall(max_pooling2d_8/PartitionedCall:output:0*
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
F__inference_flatten_2_layer_call_and_return_conditional_losses_82696442
flatten_2/PartitionedCall?
dense_4/StatefulPartitionedCallStatefulPartitionedCall"flatten_2/PartitionedCall:output:0dense_4_8269657dense_4_8269659*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_4_layer_call_and_return_conditional_losses_82696562!
dense_4/StatefulPartitionedCall?
IdentityIdentity(dense_4/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identity?
NoOpNoOp"^conv2d_14/StatefulPartitionedCall"^conv2d_15/StatefulPartitionedCall"^conv2d_16/StatefulPartitionedCall ^dense_4/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:?????????00: : : : : : : : 2F
!conv2d_14/StatefulPartitionedCall!conv2d_14/StatefulPartitionedCall2F
!conv2d_15/StatefulPartitionedCall!conv2d_15/StatefulPartitionedCall2F
!conv2d_16/StatefulPartitionedCall!conv2d_16/StatefulPartitionedCall2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall:W S
/
_output_shapes
:?????????00
 
_user_specified_nameinputs
?%
?
I__inference_sequential_4_layer_call_and_return_conditional_losses_8269890
conv2d_14_input+
conv2d_14_8269865:
conv2d_14_8269867:+
conv2d_15_8269871:
conv2d_15_8269873:+
conv2d_16_8269877:
conv2d_16_8269879:"
dense_4_8269884:	?
dense_4_8269886:
identity??!conv2d_14/StatefulPartitionedCall?!conv2d_15/StatefulPartitionedCall?!conv2d_16/StatefulPartitionedCall?dense_4/StatefulPartitionedCall?
!conv2d_14/StatefulPartitionedCallStatefulPartitionedCallconv2d_14_inputconv2d_14_8269865conv2d_14_8269867*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????00*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_14_layer_call_and_return_conditional_losses_82695802#
!conv2d_14/StatefulPartitionedCall?
max_pooling2d_6/PartitionedCallPartitionedCall*conv2d_14/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_6_layer_call_and_return_conditional_losses_82695902!
max_pooling2d_6/PartitionedCall?
!conv2d_15/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_6/PartitionedCall:output:0conv2d_15_8269871conv2d_15_8269873*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_15_layer_call_and_return_conditional_losses_82696032#
!conv2d_15/StatefulPartitionedCall?
max_pooling2d_7/PartitionedCallPartitionedCall*conv2d_15/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_7_layer_call_and_return_conditional_losses_82696132!
max_pooling2d_7/PartitionedCall?
!conv2d_16/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_7/PartitionedCall:output:0conv2d_16_8269877conv2d_16_8269879*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_16_layer_call_and_return_conditional_losses_82696262#
!conv2d_16/StatefulPartitionedCall?
max_pooling2d_8/PartitionedCallPartitionedCall*conv2d_16/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_8_layer_call_and_return_conditional_losses_82696362!
max_pooling2d_8/PartitionedCall?
flatten_2/PartitionedCallPartitionedCall(max_pooling2d_8/PartitionedCall:output:0*
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
F__inference_flatten_2_layer_call_and_return_conditional_losses_82696442
flatten_2/PartitionedCall?
dense_4/StatefulPartitionedCallStatefulPartitionedCall"flatten_2/PartitionedCall:output:0dense_4_8269884dense_4_8269886*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_4_layer_call_and_return_conditional_losses_82696562!
dense_4/StatefulPartitionedCall?
IdentityIdentity(dense_4/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identity?
NoOpNoOp"^conv2d_14/StatefulPartitionedCall"^conv2d_15/StatefulPartitionedCall"^conv2d_16/StatefulPartitionedCall ^dense_4/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:?????????00: : : : : : : : 2F
!conv2d_14/StatefulPartitionedCall!conv2d_14/StatefulPartitionedCall2F
!conv2d_15/StatefulPartitionedCall!conv2d_15/StatefulPartitionedCall2F
!conv2d_16/StatefulPartitionedCall!conv2d_16/StatefulPartitionedCall2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall:` \
/
_output_shapes
:?????????00
)
_user_specified_nameconv2d_14_input
?
?
+__inference_conv2d_15_layer_call_fn_8270087

inputs!
unknown:
	unknown_0:
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_15_layer_call_and_return_conditional_losses_82696032
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:?????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????: : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?

?
.__inference_sequential_4_layer_call_fn_8270027

inputs!
unknown:
	unknown_0:#
	unknown_1:
	unknown_2:#
	unknown_3:
	unknown_4:
	unknown_5:	?
	unknown_6:
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????**
_read_only_resource_inputs

*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_82697942
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:?????????00: : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:?????????00
 
_user_specified_nameinputs
?
h
L__inference_max_pooling2d_6_layer_call_and_return_conditional_losses_8269505

inputs
identity?
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4????????????????????????????????????*
ksize
*
paddingSAME*
strides
2	
MaxPool?
IdentityIdentityMaxPool:output:0*
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
.__inference_sequential_4_layer_call_fn_8269834
conv2d_14_input!
unknown:
	unknown_0:#
	unknown_1:
	unknown_2:#
	unknown_3:
	unknown_4:
	unknown_5:	?
	unknown_6:
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallconv2d_14_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????**
_read_only_resource_inputs

*0
config_proto 

CPU

GPU2*0J 8? *R
fMRK
I__inference_sequential_4_layer_call_and_return_conditional_losses_82697942
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:?????????00: : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:` \
/
_output_shapes
:?????????00
)
_user_specified_nameconv2d_14_input
?
M
1__inference_max_pooling2d_8_layer_call_fn_8270147

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
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_8_layer_call_and_return_conditional_losses_82696362
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?
M
1__inference_max_pooling2d_8_layer_call_fn_8270142

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
L__inference_max_pooling2d_8_layer_call_and_return_conditional_losses_82695492
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
?
G
+__inference_flatten_2_layer_call_fn_8270158

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
F__inference_flatten_2_layer_call_and_return_conditional_losses_82696442
PartitionedCallm
IdentityIdentityPartitionedCall:output:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?
M
1__inference_max_pooling2d_6_layer_call_fn_8270067

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
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_6_layer_call_and_return_conditional_losses_82695902
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????00:W S
/
_output_shapes
:?????????00
 
_user_specified_nameinputs
?
?
+__inference_conv2d_16_layer_call_fn_8270127

inputs!
unknown:
	unknown_0:
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_16_layer_call_and_return_conditional_losses_82696262
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:?????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????: : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?0
?
I__inference_sequential_4_layer_call_and_return_conditional_losses_8269949

inputsB
(conv2d_14_conv2d_readvariableop_resource:7
)conv2d_14_biasadd_readvariableop_resource:B
(conv2d_15_conv2d_readvariableop_resource:7
)conv2d_15_biasadd_readvariableop_resource:B
(conv2d_16_conv2d_readvariableop_resource:7
)conv2d_16_biasadd_readvariableop_resource:9
&dense_4_matmul_readvariableop_resource:	?5
'dense_4_biasadd_readvariableop_resource:
identity?? conv2d_14/BiasAdd/ReadVariableOp?conv2d_14/Conv2D/ReadVariableOp? conv2d_15/BiasAdd/ReadVariableOp?conv2d_15/Conv2D/ReadVariableOp? conv2d_16/BiasAdd/ReadVariableOp?conv2d_16/Conv2D/ReadVariableOp?dense_4/BiasAdd/ReadVariableOp?dense_4/MatMul/ReadVariableOp?
conv2d_14/Conv2D/ReadVariableOpReadVariableOp(conv2d_14_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02!
conv2d_14/Conv2D/ReadVariableOp?
conv2d_14/Conv2DConv2Dinputs'conv2d_14/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????00*
paddingSAME*
strides
2
conv2d_14/Conv2D?
 conv2d_14/BiasAdd/ReadVariableOpReadVariableOp)conv2d_14_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02"
 conv2d_14/BiasAdd/ReadVariableOp?
conv2d_14/BiasAddBiasAddconv2d_14/Conv2D:output:0(conv2d_14/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????002
conv2d_14/BiasAdd{
conv2d_14/EluEluconv2d_14/BiasAdd:output:0*
T0*/
_output_shapes
:?????????002
conv2d_14/Elu?
max_pooling2d_6/MaxPoolMaxPoolconv2d_14/Elu:activations:0*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2
max_pooling2d_6/MaxPool?
conv2d_15/Conv2D/ReadVariableOpReadVariableOp(conv2d_15_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02!
conv2d_15/Conv2D/ReadVariableOp?
conv2d_15/Conv2DConv2D max_pooling2d_6/MaxPool:output:0'conv2d_15/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????*
paddingSAME*
strides
2
conv2d_15/Conv2D?
 conv2d_15/BiasAdd/ReadVariableOpReadVariableOp)conv2d_15_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02"
 conv2d_15/BiasAdd/ReadVariableOp?
conv2d_15/BiasAddBiasAddconv2d_15/Conv2D:output:0(conv2d_15/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????2
conv2d_15/BiasAdd{
conv2d_15/EluEluconv2d_15/BiasAdd:output:0*
T0*/
_output_shapes
:?????????2
conv2d_15/Elu?
max_pooling2d_7/MaxPoolMaxPoolconv2d_15/Elu:activations:0*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2
max_pooling2d_7/MaxPool?
conv2d_16/Conv2D/ReadVariableOpReadVariableOp(conv2d_16_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02!
conv2d_16/Conv2D/ReadVariableOp?
conv2d_16/Conv2DConv2D max_pooling2d_7/MaxPool:output:0'conv2d_16/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????*
paddingSAME*
strides
2
conv2d_16/Conv2D?
 conv2d_16/BiasAdd/ReadVariableOpReadVariableOp)conv2d_16_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02"
 conv2d_16/BiasAdd/ReadVariableOp?
conv2d_16/BiasAddBiasAddconv2d_16/Conv2D:output:0(conv2d_16/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????2
conv2d_16/BiasAdd{
conv2d_16/EluEluconv2d_16/BiasAdd:output:0*
T0*/
_output_shapes
:?????????2
conv2d_16/Elu?
max_pooling2d_8/MaxPoolMaxPoolconv2d_16/Elu:activations:0*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2
max_pooling2d_8/MaxPools
flatten_2/ConstConst*
_output_shapes
:*
dtype0*
valueB"????   2
flatten_2/Const?
flatten_2/ReshapeReshape max_pooling2d_8/MaxPool:output:0flatten_2/Const:output:0*
T0*(
_output_shapes
:??????????2
flatten_2/Reshape?
dense_4/MatMul/ReadVariableOpReadVariableOp&dense_4_matmul_readvariableop_resource*
_output_shapes
:	?*
dtype02
dense_4/MatMul/ReadVariableOp?
dense_4/MatMulMatMulflatten_2/Reshape:output:0%dense_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
dense_4/MatMul?
dense_4/BiasAdd/ReadVariableOpReadVariableOp'dense_4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
dense_4/BiasAdd/ReadVariableOp?
dense_4/BiasAddBiasAdddense_4/MatMul:product:0&dense_4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
dense_4/BiasAdds
IdentityIdentitydense_4/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identity?
NoOpNoOp!^conv2d_14/BiasAdd/ReadVariableOp ^conv2d_14/Conv2D/ReadVariableOp!^conv2d_15/BiasAdd/ReadVariableOp ^conv2d_15/Conv2D/ReadVariableOp!^conv2d_16/BiasAdd/ReadVariableOp ^conv2d_16/Conv2D/ReadVariableOp^dense_4/BiasAdd/ReadVariableOp^dense_4/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:?????????00: : : : : : : : 2D
 conv2d_14/BiasAdd/ReadVariableOp conv2d_14/BiasAdd/ReadVariableOp2B
conv2d_14/Conv2D/ReadVariableOpconv2d_14/Conv2D/ReadVariableOp2D
 conv2d_15/BiasAdd/ReadVariableOp conv2d_15/BiasAdd/ReadVariableOp2B
conv2d_15/Conv2D/ReadVariableOpconv2d_15/Conv2D/ReadVariableOp2D
 conv2d_16/BiasAdd/ReadVariableOp conv2d_16/BiasAdd/ReadVariableOp2B
conv2d_16/Conv2D/ReadVariableOpconv2d_16/Conv2D/ReadVariableOp2@
dense_4/BiasAdd/ReadVariableOpdense_4/BiasAdd/ReadVariableOp2>
dense_4/MatMul/ReadVariableOpdense_4/MatMul/ReadVariableOp:W S
/
_output_shapes
:?????????00
 
_user_specified_nameinputs
?

?
%__inference_signature_wrapper_8269913
conv2d_14_input!
unknown:
	unknown_0:#
	unknown_1:
	unknown_2:#
	unknown_3:
	unknown_4:
	unknown_5:	?
	unknown_6:
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallconv2d_14_inputunknown	unknown_0	unknown_1	unknown_2	unknown_3	unknown_4	unknown_5	unknown_6*
Tin
2	*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????**
_read_only_resource_inputs

*0
config_proto 

CPU

GPU2*0J 8? *+
f&R$
"__inference__wrapped_model_82694962
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:?????????00: : : : : : : : 22
StatefulPartitionedCallStatefulPartitionedCall:` \
/
_output_shapes
:?????????00
)
_user_specified_nameconv2d_14_input
?
M
1__inference_max_pooling2d_7_layer_call_fn_8270107

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
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_7_layer_call_and_return_conditional_losses_82696132
PartitionedCallt
IdentityIdentityPartitionedCall:output:0*
T0*/
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?
?
F__inference_conv2d_14_layer_call_and_return_conditional_losses_8269580

inputs8
conv2d_readvariableop_resource:-
biasadd_readvariableop_resource:
identity??BiasAdd/ReadVariableOp?Conv2D/ReadVariableOp?
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
Conv2D/ReadVariableOp?
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????00*
paddingSAME*
strides
2
Conv2D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????002	
BiasAdd]
EluEluBiasAdd:output:0*
T0*/
_output_shapes
:?????????002
Elut
IdentityIdentityElu:activations:0^NoOp*
T0*/
_output_shapes
:?????????002

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????00: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:?????????00
 
_user_specified_nameinputs
?%
?
I__inference_sequential_4_layer_call_and_return_conditional_losses_8269794

inputs+
conv2d_14_8269769:
conv2d_14_8269771:+
conv2d_15_8269775:
conv2d_15_8269777:+
conv2d_16_8269781:
conv2d_16_8269783:"
dense_4_8269788:	?
dense_4_8269790:
identity??!conv2d_14/StatefulPartitionedCall?!conv2d_15/StatefulPartitionedCall?!conv2d_16/StatefulPartitionedCall?dense_4/StatefulPartitionedCall?
!conv2d_14/StatefulPartitionedCallStatefulPartitionedCallinputsconv2d_14_8269769conv2d_14_8269771*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????00*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_14_layer_call_and_return_conditional_losses_82695802#
!conv2d_14/StatefulPartitionedCall?
max_pooling2d_6/PartitionedCallPartitionedCall*conv2d_14/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_6_layer_call_and_return_conditional_losses_82695902!
max_pooling2d_6/PartitionedCall?
!conv2d_15/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_6/PartitionedCall:output:0conv2d_15_8269775conv2d_15_8269777*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_15_layer_call_and_return_conditional_losses_82696032#
!conv2d_15/StatefulPartitionedCall?
max_pooling2d_7/PartitionedCallPartitionedCall*conv2d_15/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_7_layer_call_and_return_conditional_losses_82696132!
max_pooling2d_7/PartitionedCall?
!conv2d_16/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_7/PartitionedCall:output:0conv2d_16_8269781conv2d_16_8269783*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_16_layer_call_and_return_conditional_losses_82696262#
!conv2d_16/StatefulPartitionedCall?
max_pooling2d_8/PartitionedCallPartitionedCall*conv2d_16/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_8_layer_call_and_return_conditional_losses_82696362!
max_pooling2d_8/PartitionedCall?
flatten_2/PartitionedCallPartitionedCall(max_pooling2d_8/PartitionedCall:output:0*
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
F__inference_flatten_2_layer_call_and_return_conditional_losses_82696442
flatten_2/PartitionedCall?
dense_4/StatefulPartitionedCallStatefulPartitionedCall"flatten_2/PartitionedCall:output:0dense_4_8269788dense_4_8269790*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_4_layer_call_and_return_conditional_losses_82696562!
dense_4/StatefulPartitionedCall?
IdentityIdentity(dense_4/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identity?
NoOpNoOp"^conv2d_14/StatefulPartitionedCall"^conv2d_15/StatefulPartitionedCall"^conv2d_16/StatefulPartitionedCall ^dense_4/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:?????????00: : : : : : : : 2F
!conv2d_14/StatefulPartitionedCall!conv2d_14/StatefulPartitionedCall2F
!conv2d_15/StatefulPartitionedCall!conv2d_15/StatefulPartitionedCall2F
!conv2d_16/StatefulPartitionedCall!conv2d_16/StatefulPartitionedCall2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall:W S
/
_output_shapes
:?????????00
 
_user_specified_nameinputs
?
h
L__inference_max_pooling2d_8_layer_call_and_return_conditional_losses_8270137

inputs
identity?
MaxPoolMaxPoolinputs*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2	
MaxPooll
IdentityIdentityMaxPool:output:0*
T0*/
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?

?
D__inference_dense_4_layer_call_and_return_conditional_losses_8269656

inputs1
matmul_readvariableop_resource:	?-
biasadd_readvariableop_resource:
identity??BiasAdd/ReadVariableOp?MatMul/ReadVariableOp?
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	?*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2	
BiasAddk
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:??????????: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
h
L__inference_max_pooling2d_7_layer_call_and_return_conditional_losses_8270097

inputs
identity?
MaxPoolMaxPoolinputs*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2	
MaxPooll
IdentityIdentityMaxPool:output:0*
T0*/
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?
?
)__inference_dense_4_layer_call_fn_8270177

inputs
unknown:	?
	unknown_0:
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_4_layer_call_and_return_conditional_losses_82696562
StatefulPartitionedCall{
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:??????????: : 22
StatefulPartitionedCallStatefulPartitionedCall:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?

?
D__inference_dense_4_layer_call_and_return_conditional_losses_8270168

inputs1
matmul_readvariableop_resource:	?-
biasadd_readvariableop_resource:
identity??BiasAdd/ReadVariableOp?MatMul/ReadVariableOp?
MatMul/ReadVariableOpReadVariableOpmatmul_readvariableop_resource*
_output_shapes
:	?*
dtype02
MatMul/ReadVariableOps
MatMulMatMulinputsMatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
MatMul?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddMatMul:product:0BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2	
BiasAddk
IdentityIdentityBiasAdd:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*+
_input_shapes
:??????????: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
MatMul/ReadVariableOpMatMul/ReadVariableOp:P L
(
_output_shapes
:??????????
 
_user_specified_nameinputs
?
M
1__inference_max_pooling2d_6_layer_call_fn_8270062

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
L__inference_max_pooling2d_6_layer_call_and_return_conditional_losses_82695052
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
?
h
L__inference_max_pooling2d_8_layer_call_and_return_conditional_losses_8270132

inputs
identity?
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4????????????????????????????????????*
ksize
*
paddingSAME*
strides
2	
MaxPool?
IdentityIdentityMaxPool:output:0*
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
?
b
F__inference_flatten_2_layer_call_and_return_conditional_losses_8270153

inputs
identity_
ConstConst*
_output_shapes
:*
dtype0*
valueB"????   2
Consth
ReshapeReshapeinputsConst:output:0*
T0*(
_output_shapes
:??????????2	
Reshapee
IdentityIdentityReshape:output:0*
T0*(
_output_shapes
:??????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?'
?
#__inference__traced_restore_8270258
file_prefix;
!assignvariableop_conv2d_14_kernel:/
!assignvariableop_1_conv2d_14_bias:=
#assignvariableop_2_conv2d_15_kernel:/
!assignvariableop_3_conv2d_15_bias:=
#assignvariableop_4_conv2d_16_kernel:/
!assignvariableop_5_conv2d_16_bias:4
!assignvariableop_6_dense_4_kernel:	?-
assignvariableop_7_dense_4_bias:

identity_9??AssignVariableOp?AssignVariableOp_1?AssignVariableOp_2?AssignVariableOp_3?AssignVariableOp_4?AssignVariableOp_5?AssignVariableOp_6?AssignVariableOp_7?
RestoreV2/tensor_namesConst"/device:CPU:0*
_output_shapes
:	*
dtype0*?
value?B?	B6layer_with_weights-0/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-0/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-1/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-1/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-2/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-2/bias/.ATTRIBUTES/VARIABLE_VALUEB6layer_with_weights-3/kernel/.ATTRIBUTES/VARIABLE_VALUEB4layer_with_weights-3/bias/.ATTRIBUTES/VARIABLE_VALUEB_CHECKPOINTABLE_OBJECT_GRAPH2
RestoreV2/tensor_names?
RestoreV2/shape_and_slicesConst"/device:CPU:0*
_output_shapes
:	*
dtype0*%
valueB	B B B B B B B B B 2
RestoreV2/shape_and_slices?
	RestoreV2	RestoreV2file_prefixRestoreV2/tensor_names:output:0#RestoreV2/shape_and_slices:output:0"/device:CPU:0*8
_output_shapes&
$:::::::::*
dtypes
2	2
	RestoreV2g
IdentityIdentityRestoreV2:tensors:0"/device:CPU:0*
T0*
_output_shapes
:2

Identity?
AssignVariableOpAssignVariableOp!assignvariableop_conv2d_14_kernelIdentity:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOpk

Identity_1IdentityRestoreV2:tensors:1"/device:CPU:0*
T0*
_output_shapes
:2

Identity_1?
AssignVariableOp_1AssignVariableOp!assignvariableop_1_conv2d_14_biasIdentity_1:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_1k

Identity_2IdentityRestoreV2:tensors:2"/device:CPU:0*
T0*
_output_shapes
:2

Identity_2?
AssignVariableOp_2AssignVariableOp#assignvariableop_2_conv2d_15_kernelIdentity_2:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_2k

Identity_3IdentityRestoreV2:tensors:3"/device:CPU:0*
T0*
_output_shapes
:2

Identity_3?
AssignVariableOp_3AssignVariableOp!assignvariableop_3_conv2d_15_biasIdentity_3:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_3k

Identity_4IdentityRestoreV2:tensors:4"/device:CPU:0*
T0*
_output_shapes
:2

Identity_4?
AssignVariableOp_4AssignVariableOp#assignvariableop_4_conv2d_16_kernelIdentity_4:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_4k

Identity_5IdentityRestoreV2:tensors:5"/device:CPU:0*
T0*
_output_shapes
:2

Identity_5?
AssignVariableOp_5AssignVariableOp!assignvariableop_5_conv2d_16_biasIdentity_5:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_5k

Identity_6IdentityRestoreV2:tensors:6"/device:CPU:0*
T0*
_output_shapes
:2

Identity_6?
AssignVariableOp_6AssignVariableOp!assignvariableop_6_dense_4_kernelIdentity_6:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_6k

Identity_7IdentityRestoreV2:tensors:7"/device:CPU:0*
T0*
_output_shapes
:2

Identity_7?
AssignVariableOp_7AssignVariableOpassignvariableop_7_dense_4_biasIdentity_7:output:0"/device:CPU:0*
_output_shapes
 *
dtype02
AssignVariableOp_79
NoOpNoOp"/device:CPU:0*
_output_shapes
 2
NoOp?

Identity_8Identityfile_prefix^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_2^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7^NoOp"/device:CPU:0*
T0*
_output_shapes
: 2

Identity_8c

Identity_9IdentityIdentity_8:output:0^NoOp_1*
T0*
_output_shapes
: 2

Identity_9?
NoOp_1NoOp^AssignVariableOp^AssignVariableOp_1^AssignVariableOp_2^AssignVariableOp_3^AssignVariableOp_4^AssignVariableOp_5^AssignVariableOp_6^AssignVariableOp_7*"
_acd_function_control_output(*
_output_shapes
 2
NoOp_1"!

identity_9Identity_9:output:0*%
_input_shapes
: : : : : : : : : 2$
AssignVariableOpAssignVariableOp2(
AssignVariableOp_1AssignVariableOp_12(
AssignVariableOp_2AssignVariableOp_22(
AssignVariableOp_3AssignVariableOp_32(
AssignVariableOp_4AssignVariableOp_42(
AssignVariableOp_5AssignVariableOp_52(
AssignVariableOp_6AssignVariableOp_62(
AssignVariableOp_7AssignVariableOp_7:C ?

_output_shapes
: 
%
_user_specified_namefile_prefix
?
h
L__inference_max_pooling2d_8_layer_call_and_return_conditional_losses_8269636

inputs
identity?
MaxPoolMaxPoolinputs*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2	
MaxPooll
IdentityIdentityMaxPool:output:0*
T0*/
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?
?
F__inference_conv2d_16_layer_call_and_return_conditional_losses_8269626

inputs8
conv2d_readvariableop_resource:-
biasadd_readvariableop_resource:
identity??BiasAdd/ReadVariableOp?Conv2D/ReadVariableOp?
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
Conv2D/ReadVariableOp?
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????*
paddingSAME*
strides
2
Conv2D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????2	
BiasAdd]
EluEluBiasAdd:output:0*
T0*/
_output_shapes
:?????????2
Elut
IdentityIdentityElu:activations:0^NoOp*
T0*/
_output_shapes
:?????????2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?
h
L__inference_max_pooling2d_6_layer_call_and_return_conditional_losses_8269590

inputs
identity?
MaxPoolMaxPoolinputs*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2	
MaxPooll
IdentityIdentityMaxPool:output:0*
T0*/
_output_shapes
:?????????2

Identity"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*.
_input_shapes
:?????????00:W S
/
_output_shapes
:?????????00
 
_user_specified_nameinputs
?
?
F__inference_conv2d_15_layer_call_and_return_conditional_losses_8269603

inputs8
conv2d_readvariableop_resource:-
biasadd_readvariableop_resource:
identity??BiasAdd/ReadVariableOp?Conv2D/ReadVariableOp?
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
Conv2D/ReadVariableOp?
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????*
paddingSAME*
strides
2
Conv2D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????2	
BiasAdd]
EluEluBiasAdd:output:0*
T0*/
_output_shapes
:?????????2
Elut
IdentityIdentityElu:activations:0^NoOp*
T0*/
_output_shapes
:?????????2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?%
?
I__inference_sequential_4_layer_call_and_return_conditional_losses_8269862
conv2d_14_input+
conv2d_14_8269837:
conv2d_14_8269839:+
conv2d_15_8269843:
conv2d_15_8269845:+
conv2d_16_8269849:
conv2d_16_8269851:"
dense_4_8269856:	?
dense_4_8269858:
identity??!conv2d_14/StatefulPartitionedCall?!conv2d_15/StatefulPartitionedCall?!conv2d_16/StatefulPartitionedCall?dense_4/StatefulPartitionedCall?
!conv2d_14/StatefulPartitionedCallStatefulPartitionedCallconv2d_14_inputconv2d_14_8269837conv2d_14_8269839*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????00*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_14_layer_call_and_return_conditional_losses_82695802#
!conv2d_14/StatefulPartitionedCall?
max_pooling2d_6/PartitionedCallPartitionedCall*conv2d_14/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_6_layer_call_and_return_conditional_losses_82695902!
max_pooling2d_6/PartitionedCall?
!conv2d_15/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_6/PartitionedCall:output:0conv2d_15_8269843conv2d_15_8269845*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_15_layer_call_and_return_conditional_losses_82696032#
!conv2d_15/StatefulPartitionedCall?
max_pooling2d_7/PartitionedCallPartitionedCall*conv2d_15/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_7_layer_call_and_return_conditional_losses_82696132!
max_pooling2d_7/PartitionedCall?
!conv2d_16/StatefulPartitionedCallStatefulPartitionedCall(max_pooling2d_7/PartitionedCall:output:0conv2d_16_8269849conv2d_16_8269851*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_16_layer_call_and_return_conditional_losses_82696262#
!conv2d_16/StatefulPartitionedCall?
max_pooling2d_8/PartitionedCallPartitionedCall*conv2d_16/StatefulPartitionedCall:output:0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????* 
_read_only_resource_inputs
 *0
config_proto 

CPU

GPU2*0J 8? *U
fPRN
L__inference_max_pooling2d_8_layer_call_and_return_conditional_losses_82696362!
max_pooling2d_8/PartitionedCall?
flatten_2/PartitionedCallPartitionedCall(max_pooling2d_8/PartitionedCall:output:0*
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
F__inference_flatten_2_layer_call_and_return_conditional_losses_82696442
flatten_2/PartitionedCall?
dense_4/StatefulPartitionedCallStatefulPartitionedCall"flatten_2/PartitionedCall:output:0dense_4_8269856dense_4_8269858*
Tin
2*
Tout
2*
_collective_manager_ids
 *'
_output_shapes
:?????????*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *M
fHRF
D__inference_dense_4_layer_call_and_return_conditional_losses_82696562!
dense_4/StatefulPartitionedCall?
IdentityIdentity(dense_4/StatefulPartitionedCall:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identity?
NoOpNoOp"^conv2d_14/StatefulPartitionedCall"^conv2d_15/StatefulPartitionedCall"^conv2d_16/StatefulPartitionedCall ^dense_4/StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:?????????00: : : : : : : : 2F
!conv2d_14/StatefulPartitionedCall!conv2d_14/StatefulPartitionedCall2F
!conv2d_15/StatefulPartitionedCall!conv2d_15/StatefulPartitionedCall2F
!conv2d_16/StatefulPartitionedCall!conv2d_16/StatefulPartitionedCall2B
dense_4/StatefulPartitionedCalldense_4/StatefulPartitionedCall:` \
/
_output_shapes
:?????????00
)
_user_specified_nameconv2d_14_input
?
h
L__inference_max_pooling2d_6_layer_call_and_return_conditional_losses_8270052

inputs
identity?
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4????????????????????????????????????*
ksize
*
paddingSAME*
strides
2	
MaxPool?
IdentityIdentityMaxPool:output:0*
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
F__inference_conv2d_15_layer_call_and_return_conditional_losses_8270078

inputs8
conv2d_readvariableop_resource:-
biasadd_readvariableop_resource:
identity??BiasAdd/ReadVariableOp?Conv2D/ReadVariableOp?
Conv2D/ReadVariableOpReadVariableOpconv2d_readvariableop_resource*&
_output_shapes
:*
dtype02
Conv2D/ReadVariableOp?
Conv2DConv2DinputsConv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????*
paddingSAME*
strides
2
Conv2D?
BiasAdd/ReadVariableOpReadVariableOpbiasadd_readvariableop_resource*
_output_shapes
:*
dtype02
BiasAdd/ReadVariableOp?
BiasAddBiasAddConv2D:output:0BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????2	
BiasAdd]
EluEluBiasAdd:output:0*
T0*/
_output_shapes
:?????????2
Elut
IdentityIdentityElu:activations:0^NoOp*
T0*/
_output_shapes
:?????????2

Identity
NoOpNoOp^BiasAdd/ReadVariableOp^Conv2D/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????: : 20
BiasAdd/ReadVariableOpBiasAdd/ReadVariableOp2.
Conv2D/ReadVariableOpConv2D/ReadVariableOp:W S
/
_output_shapes
:?????????
 
_user_specified_nameinputs
?
h
L__inference_max_pooling2d_7_layer_call_and_return_conditional_losses_8270092

inputs
identity?
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4????????????????????????????????????*
ksize
*
paddingSAME*
strides
2	
MaxPool?
IdentityIdentityMaxPool:output:0*
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
?
h
L__inference_max_pooling2d_8_layer_call_and_return_conditional_losses_8269549

inputs
identity?
MaxPoolMaxPoolinputs*J
_output_shapes8
6:4????????????????????????????????????*
ksize
*
paddingSAME*
strides
2	
MaxPool?
IdentityIdentityMaxPool:output:0*
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
?0
?
I__inference_sequential_4_layer_call_and_return_conditional_losses_8269985

inputsB
(conv2d_14_conv2d_readvariableop_resource:7
)conv2d_14_biasadd_readvariableop_resource:B
(conv2d_15_conv2d_readvariableop_resource:7
)conv2d_15_biasadd_readvariableop_resource:B
(conv2d_16_conv2d_readvariableop_resource:7
)conv2d_16_biasadd_readvariableop_resource:9
&dense_4_matmul_readvariableop_resource:	?5
'dense_4_biasadd_readvariableop_resource:
identity?? conv2d_14/BiasAdd/ReadVariableOp?conv2d_14/Conv2D/ReadVariableOp? conv2d_15/BiasAdd/ReadVariableOp?conv2d_15/Conv2D/ReadVariableOp? conv2d_16/BiasAdd/ReadVariableOp?conv2d_16/Conv2D/ReadVariableOp?dense_4/BiasAdd/ReadVariableOp?dense_4/MatMul/ReadVariableOp?
conv2d_14/Conv2D/ReadVariableOpReadVariableOp(conv2d_14_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02!
conv2d_14/Conv2D/ReadVariableOp?
conv2d_14/Conv2DConv2Dinputs'conv2d_14/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????00*
paddingSAME*
strides
2
conv2d_14/Conv2D?
 conv2d_14/BiasAdd/ReadVariableOpReadVariableOp)conv2d_14_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02"
 conv2d_14/BiasAdd/ReadVariableOp?
conv2d_14/BiasAddBiasAddconv2d_14/Conv2D:output:0(conv2d_14/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????002
conv2d_14/BiasAdd{
conv2d_14/EluEluconv2d_14/BiasAdd:output:0*
T0*/
_output_shapes
:?????????002
conv2d_14/Elu?
max_pooling2d_6/MaxPoolMaxPoolconv2d_14/Elu:activations:0*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2
max_pooling2d_6/MaxPool?
conv2d_15/Conv2D/ReadVariableOpReadVariableOp(conv2d_15_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02!
conv2d_15/Conv2D/ReadVariableOp?
conv2d_15/Conv2DConv2D max_pooling2d_6/MaxPool:output:0'conv2d_15/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????*
paddingSAME*
strides
2
conv2d_15/Conv2D?
 conv2d_15/BiasAdd/ReadVariableOpReadVariableOp)conv2d_15_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02"
 conv2d_15/BiasAdd/ReadVariableOp?
conv2d_15/BiasAddBiasAddconv2d_15/Conv2D:output:0(conv2d_15/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????2
conv2d_15/BiasAdd{
conv2d_15/EluEluconv2d_15/BiasAdd:output:0*
T0*/
_output_shapes
:?????????2
conv2d_15/Elu?
max_pooling2d_7/MaxPoolMaxPoolconv2d_15/Elu:activations:0*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2
max_pooling2d_7/MaxPool?
conv2d_16/Conv2D/ReadVariableOpReadVariableOp(conv2d_16_conv2d_readvariableop_resource*&
_output_shapes
:*
dtype02!
conv2d_16/Conv2D/ReadVariableOp?
conv2d_16/Conv2DConv2D max_pooling2d_7/MaxPool:output:0'conv2d_16/Conv2D/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????*
paddingSAME*
strides
2
conv2d_16/Conv2D?
 conv2d_16/BiasAdd/ReadVariableOpReadVariableOp)conv2d_16_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02"
 conv2d_16/BiasAdd/ReadVariableOp?
conv2d_16/BiasAddBiasAddconv2d_16/Conv2D:output:0(conv2d_16/BiasAdd/ReadVariableOp:value:0*
T0*/
_output_shapes
:?????????2
conv2d_16/BiasAdd{
conv2d_16/EluEluconv2d_16/BiasAdd:output:0*
T0*/
_output_shapes
:?????????2
conv2d_16/Elu?
max_pooling2d_8/MaxPoolMaxPoolconv2d_16/Elu:activations:0*/
_output_shapes
:?????????*
ksize
*
paddingSAME*
strides
2
max_pooling2d_8/MaxPools
flatten_2/ConstConst*
_output_shapes
:*
dtype0*
valueB"????   2
flatten_2/Const?
flatten_2/ReshapeReshape max_pooling2d_8/MaxPool:output:0flatten_2/Const:output:0*
T0*(
_output_shapes
:??????????2
flatten_2/Reshape?
dense_4/MatMul/ReadVariableOpReadVariableOp&dense_4_matmul_readvariableop_resource*
_output_shapes
:	?*
dtype02
dense_4/MatMul/ReadVariableOp?
dense_4/MatMulMatMulflatten_2/Reshape:output:0%dense_4/MatMul/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
dense_4/MatMul?
dense_4/BiasAdd/ReadVariableOpReadVariableOp'dense_4_biasadd_readvariableop_resource*
_output_shapes
:*
dtype02 
dense_4/BiasAdd/ReadVariableOp?
dense_4/BiasAddBiasAdddense_4/MatMul:product:0&dense_4/BiasAdd/ReadVariableOp:value:0*
T0*'
_output_shapes
:?????????2
dense_4/BiasAdds
IdentityIdentitydense_4/BiasAdd:output:0^NoOp*
T0*'
_output_shapes
:?????????2

Identity?
NoOpNoOp!^conv2d_14/BiasAdd/ReadVariableOp ^conv2d_14/Conv2D/ReadVariableOp!^conv2d_15/BiasAdd/ReadVariableOp ^conv2d_15/Conv2D/ReadVariableOp!^conv2d_16/BiasAdd/ReadVariableOp ^conv2d_16/Conv2D/ReadVariableOp^dense_4/BiasAdd/ReadVariableOp^dense_4/MatMul/ReadVariableOp*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*>
_input_shapes-
+:?????????00: : : : : : : : 2D
 conv2d_14/BiasAdd/ReadVariableOp conv2d_14/BiasAdd/ReadVariableOp2B
conv2d_14/Conv2D/ReadVariableOpconv2d_14/Conv2D/ReadVariableOp2D
 conv2d_15/BiasAdd/ReadVariableOp conv2d_15/BiasAdd/ReadVariableOp2B
conv2d_15/Conv2D/ReadVariableOpconv2d_15/Conv2D/ReadVariableOp2D
 conv2d_16/BiasAdd/ReadVariableOp conv2d_16/BiasAdd/ReadVariableOp2B
conv2d_16/Conv2D/ReadVariableOpconv2d_16/Conv2D/ReadVariableOp2@
dense_4/BiasAdd/ReadVariableOpdense_4/BiasAdd/ReadVariableOp2>
dense_4/MatMul/ReadVariableOpdense_4/MatMul/ReadVariableOp:W S
/
_output_shapes
:?????????00
 
_user_specified_nameinputs
?
M
1__inference_max_pooling2d_7_layer_call_fn_8270102

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
L__inference_max_pooling2d_7_layer_call_and_return_conditional_losses_82695272
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
?
?
+__inference_conv2d_14_layer_call_fn_8270047

inputs!
unknown:
	unknown_0:
identity??StatefulPartitionedCall?
StatefulPartitionedCallStatefulPartitionedCallinputsunknown	unknown_0*
Tin
2*
Tout
2*
_collective_manager_ids
 */
_output_shapes
:?????????00*$
_read_only_resource_inputs
*0
config_proto 

CPU

GPU2*0J 8? *O
fJRH
F__inference_conv2d_14_layer_call_and_return_conditional_losses_82695802
StatefulPartitionedCall?
IdentityIdentity StatefulPartitionedCall:output:0^NoOp*
T0*/
_output_shapes
:?????????002

Identityh
NoOpNoOp^StatefulPartitionedCall*"
_acd_function_control_output(*
_output_shapes
 2
NoOp"
identityIdentity:output:0*(
_construction_contextkEagerRuntime*2
_input_shapes!
:?????????00: : 22
StatefulPartitionedCallStatefulPartitionedCall:W S
/
_output_shapes
:?????????00
 
_user_specified_nameinputs"?L
saver_filename:0StatefulPartitionedCall_1:0StatefulPartitionedCall_28"
saved_model_main_op

NoOp*>
__saved_model_init_op%#
__saved_model_init_op

NoOp*?
serving_default?
S
conv2d_14_input@
!serving_default_conv2d_14_input:0?????????00;
dense_40
StatefulPartitionedCall:0?????????tensorflow/serving/predict:??
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
layer-6
layer_with_weights-3
layer-7
	regularization_losses

trainable_variables
	variables
	keras_api

signatures
c_default_save_signature
*d&call_and_return_all_conditional_losses
e__call__"
_tf_keras_sequential
?

kernel
bias
regularization_losses
trainable_variables
	variables
	keras_api
*f&call_and_return_all_conditional_losses
g__call__"
_tf_keras_layer
?
regularization_losses
trainable_variables
	variables
	keras_api
*h&call_and_return_all_conditional_losses
i__call__"
_tf_keras_layer
?

kernel
bias
regularization_losses
trainable_variables
	variables
	keras_api
*j&call_and_return_all_conditional_losses
k__call__"
_tf_keras_layer
?
regularization_losses
trainable_variables
 	variables
!	keras_api
*l&call_and_return_all_conditional_losses
m__call__"
_tf_keras_layer
?

"kernel
#bias
$regularization_losses
%trainable_variables
&	variables
'	keras_api
*n&call_and_return_all_conditional_losses
o__call__"
_tf_keras_layer
?
(regularization_losses
)trainable_variables
*	variables
+	keras_api
*p&call_and_return_all_conditional_losses
q__call__"
_tf_keras_layer
?
,regularization_losses
-trainable_variables
.	variables
/	keras_api
*r&call_and_return_all_conditional_losses
s__call__"
_tf_keras_layer
?

0kernel
1bias
2regularization_losses
3trainable_variables
4	variables
5	keras_api
*t&call_and_return_all_conditional_losses
u__call__"
_tf_keras_layer
 "
trackable_list_wrapper
X
0
1
2
3
"4
#5
06
17"
trackable_list_wrapper
X
0
1
2
3
"4
#5
06
17"
trackable_list_wrapper
?
	regularization_losses
6layer_regularization_losses

trainable_variables
	variables

7layers
8non_trainable_variables
9layer_metrics
:metrics
e__call__
c_default_save_signature
*d&call_and_return_all_conditional_losses
&d"call_and_return_conditional_losses"
_generic_user_object
,
vserving_default"
signature_map
*:(2conv2d_14/kernel
:2conv2d_14/bias
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
?
regularization_losses
;layer_regularization_losses
trainable_variables
	variables

<layers
=non_trainable_variables
>layer_metrics
?metrics
g__call__
*f&call_and_return_all_conditional_losses
&f"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
regularization_losses
@layer_regularization_losses
trainable_variables
	variables

Alayers
Bnon_trainable_variables
Clayer_metrics
Dmetrics
i__call__
*h&call_and_return_all_conditional_losses
&h"call_and_return_conditional_losses"
_generic_user_object
*:(2conv2d_15/kernel
:2conv2d_15/bias
 "
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
.
0
1"
trackable_list_wrapper
?
regularization_losses
Elayer_regularization_losses
trainable_variables
	variables

Flayers
Gnon_trainable_variables
Hlayer_metrics
Imetrics
k__call__
*j&call_and_return_all_conditional_losses
&j"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
regularization_losses
Jlayer_regularization_losses
trainable_variables
 	variables

Klayers
Lnon_trainable_variables
Mlayer_metrics
Nmetrics
m__call__
*l&call_and_return_all_conditional_losses
&l"call_and_return_conditional_losses"
_generic_user_object
*:(2conv2d_16/kernel
:2conv2d_16/bias
 "
trackable_list_wrapper
.
"0
#1"
trackable_list_wrapper
.
"0
#1"
trackable_list_wrapper
?
$regularization_losses
Olayer_regularization_losses
%trainable_variables
&	variables

Players
Qnon_trainable_variables
Rlayer_metrics
Smetrics
o__call__
*n&call_and_return_all_conditional_losses
&n"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
(regularization_losses
Tlayer_regularization_losses
)trainable_variables
*	variables

Ulayers
Vnon_trainable_variables
Wlayer_metrics
Xmetrics
q__call__
*p&call_and_return_all_conditional_losses
&p"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
 "
trackable_list_wrapper
 "
trackable_list_wrapper
?
,regularization_losses
Ylayer_regularization_losses
-trainable_variables
.	variables

Zlayers
[non_trainable_variables
\layer_metrics
]metrics
s__call__
*r&call_and_return_all_conditional_losses
&r"call_and_return_conditional_losses"
_generic_user_object
!:	?2dense_4/kernel
:2dense_4/bias
 "
trackable_list_wrapper
.
00
11"
trackable_list_wrapper
.
00
11"
trackable_list_wrapper
?
2regularization_losses
^layer_regularization_losses
3trainable_variables
4	variables

_layers
`non_trainable_variables
alayer_metrics
bmetrics
u__call__
*t&call_and_return_all_conditional_losses
&t"call_and_return_conditional_losses"
_generic_user_object
 "
trackable_list_wrapper
X
0
1
2
3
4
5
6
7"
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
?B?
"__inference__wrapped_model_8269496conv2d_14_input"?
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
I__inference_sequential_4_layer_call_and_return_conditional_losses_8269949
I__inference_sequential_4_layer_call_and_return_conditional_losses_8269985
I__inference_sequential_4_layer_call_and_return_conditional_losses_8269862
I__inference_sequential_4_layer_call_and_return_conditional_losses_8269890?
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
?2?
.__inference_sequential_4_layer_call_fn_8269682
.__inference_sequential_4_layer_call_fn_8270006
.__inference_sequential_4_layer_call_fn_8270027
.__inference_sequential_4_layer_call_fn_8269834?
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
F__inference_conv2d_14_layer_call_and_return_conditional_losses_8270038?
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
+__inference_conv2d_14_layer_call_fn_8270047?
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
L__inference_max_pooling2d_6_layer_call_and_return_conditional_losses_8270052
L__inference_max_pooling2d_6_layer_call_and_return_conditional_losses_8270057?
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
1__inference_max_pooling2d_6_layer_call_fn_8270062
1__inference_max_pooling2d_6_layer_call_fn_8270067?
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
F__inference_conv2d_15_layer_call_and_return_conditional_losses_8270078?
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
+__inference_conv2d_15_layer_call_fn_8270087?
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
L__inference_max_pooling2d_7_layer_call_and_return_conditional_losses_8270092
L__inference_max_pooling2d_7_layer_call_and_return_conditional_losses_8270097?
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
1__inference_max_pooling2d_7_layer_call_fn_8270102
1__inference_max_pooling2d_7_layer_call_fn_8270107?
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
F__inference_conv2d_16_layer_call_and_return_conditional_losses_8270118?
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
+__inference_conv2d_16_layer_call_fn_8270127?
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
L__inference_max_pooling2d_8_layer_call_and_return_conditional_losses_8270132
L__inference_max_pooling2d_8_layer_call_and_return_conditional_losses_8270137?
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
1__inference_max_pooling2d_8_layer_call_fn_8270142
1__inference_max_pooling2d_8_layer_call_fn_8270147?
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
F__inference_flatten_2_layer_call_and_return_conditional_losses_8270153?
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
+__inference_flatten_2_layer_call_fn_8270158?
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
D__inference_dense_4_layer_call_and_return_conditional_losses_8270168?
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
)__inference_dense_4_layer_call_fn_8270177?
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
?B?
%__inference_signature_wrapper_8269913conv2d_14_input"?
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
"__inference__wrapped_model_8269496"#01@?=
6?3
1?.
conv2d_14_input?????????00
? "1?.
,
dense_4!?
dense_4??????????
F__inference_conv2d_14_layer_call_and_return_conditional_losses_8270038l7?4
-?*
(?%
inputs?????????00
? "-?*
#? 
0?????????00
? ?
+__inference_conv2d_14_layer_call_fn_8270047_7?4
-?*
(?%
inputs?????????00
? " ??????????00?
F__inference_conv2d_15_layer_call_and_return_conditional_losses_8270078l7?4
-?*
(?%
inputs?????????
? "-?*
#? 
0?????????
? ?
+__inference_conv2d_15_layer_call_fn_8270087_7?4
-?*
(?%
inputs?????????
? " ???????????
F__inference_conv2d_16_layer_call_and_return_conditional_losses_8270118l"#7?4
-?*
(?%
inputs?????????
? "-?*
#? 
0?????????
? ?
+__inference_conv2d_16_layer_call_fn_8270127_"#7?4
-?*
(?%
inputs?????????
? " ???????????
D__inference_dense_4_layer_call_and_return_conditional_losses_8270168]010?-
&?#
!?
inputs??????????
? "%?"
?
0?????????
? }
)__inference_dense_4_layer_call_fn_8270177P010?-
&?#
!?
inputs??????????
? "???????????
F__inference_flatten_2_layer_call_and_return_conditional_losses_8270153a7?4
-?*
(?%
inputs?????????
? "&?#
?
0??????????
? ?
+__inference_flatten_2_layer_call_fn_8270158T7?4
-?*
(?%
inputs?????????
? "????????????
L__inference_max_pooling2d_6_layer_call_and_return_conditional_losses_8270052?R?O
H?E
C?@
inputs4????????????????????????????????????
? "H?E
>?;
04????????????????????????????????????
? ?
L__inference_max_pooling2d_6_layer_call_and_return_conditional_losses_8270057h7?4
-?*
(?%
inputs?????????00
? "-?*
#? 
0?????????
? ?
1__inference_max_pooling2d_6_layer_call_fn_8270062?R?O
H?E
C?@
inputs4????????????????????????????????????
? ";?84?????????????????????????????????????
1__inference_max_pooling2d_6_layer_call_fn_8270067[7?4
-?*
(?%
inputs?????????00
? " ???????????
L__inference_max_pooling2d_7_layer_call_and_return_conditional_losses_8270092?R?O
H?E
C?@
inputs4????????????????????????????????????
? "H?E
>?;
04????????????????????????????????????
? ?
L__inference_max_pooling2d_7_layer_call_and_return_conditional_losses_8270097h7?4
-?*
(?%
inputs?????????
? "-?*
#? 
0?????????
? ?
1__inference_max_pooling2d_7_layer_call_fn_8270102?R?O
H?E
C?@
inputs4????????????????????????????????????
? ";?84?????????????????????????????????????
1__inference_max_pooling2d_7_layer_call_fn_8270107[7?4
-?*
(?%
inputs?????????
? " ???????????
L__inference_max_pooling2d_8_layer_call_and_return_conditional_losses_8270132?R?O
H?E
C?@
inputs4????????????????????????????????????
? "H?E
>?;
04????????????????????????????????????
? ?
L__inference_max_pooling2d_8_layer_call_and_return_conditional_losses_8270137h7?4
-?*
(?%
inputs?????????
? "-?*
#? 
0?????????
? ?
1__inference_max_pooling2d_8_layer_call_fn_8270142?R?O
H?E
C?@
inputs4????????????????????????????????????
? ";?84?????????????????????????????????????
1__inference_max_pooling2d_8_layer_call_fn_8270147[7?4
-?*
(?%
inputs?????????
? " ???????????
I__inference_sequential_4_layer_call_and_return_conditional_losses_8269862{"#01H?E
>?;
1?.
conv2d_14_input?????????00
p 

 
? "%?"
?
0?????????
? ?
I__inference_sequential_4_layer_call_and_return_conditional_losses_8269890{"#01H?E
>?;
1?.
conv2d_14_input?????????00
p

 
? "%?"
?
0?????????
? ?
I__inference_sequential_4_layer_call_and_return_conditional_losses_8269949r"#01??<
5?2
(?%
inputs?????????00
p 

 
? "%?"
?
0?????????
? ?
I__inference_sequential_4_layer_call_and_return_conditional_losses_8269985r"#01??<
5?2
(?%
inputs?????????00
p

 
? "%?"
?
0?????????
? ?
.__inference_sequential_4_layer_call_fn_8269682n"#01H?E
>?;
1?.
conv2d_14_input?????????00
p 

 
? "???????????
.__inference_sequential_4_layer_call_fn_8269834n"#01H?E
>?;
1?.
conv2d_14_input?????????00
p

 
? "???????????
.__inference_sequential_4_layer_call_fn_8270006e"#01??<
5?2
(?%
inputs?????????00
p 

 
? "???????????
.__inference_sequential_4_layer_call_fn_8270027e"#01??<
5?2
(?%
inputs?????????00
p

 
? "???????????
%__inference_signature_wrapper_8269913?"#01S?P
? 
I?F
D
conv2d_14_input1?.
conv2d_14_input?????????00"1?.
,
dense_4!?
dense_4?????????