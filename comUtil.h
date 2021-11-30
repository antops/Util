/*
*	created by liyun 2018/5/24
*   function ����
*   version 1.0
*/


#ifndef COMUTIL_H
#define COMUTIL_H

#include "GraphTrans.h"
#include <vtkjsoncpp/json/json.h>

void packGraphTransToJson(const GraphTrans & graphTrans, Json::Value & js);

bool parseJsonToGraphTrans(const Json::Value & js, GraphTrans & graphTrans, double scale = 1);


#endif // COMUTIL_H
