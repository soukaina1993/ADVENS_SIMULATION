//
// Created by cornelia.blanke on 10.06.2025.
//

#ifndef FLUIDS_DBCONNECTION_H
#define FLUIDS_DBCONNECTION_H

#include<string>
#include<vector>
#include <pqxx/pqxx>

using namespace std;


class DBConnection {
public:
    DBConnection(const char* dbname, const char* user, const char* password,
                 const char* host="localhost", int port=5432);
    ~DBConnection() { delete conn; };

private:
    pqxx::connection* conn=nullptr;
    string dbname, user, password, host="localhost";        // database
    int port=5432;

public:
    pqxx::connection* Conn() const { return conn; }
    string query_DB() const;
};


class DBSchema {
public:
    DBSchema(const char* schemaname);     // if true, everything in the schema will be deleted
    ~DBSchema()=default;

private:
    string name;

public:
    string SchemaName() const { return name; }

    string query_create_schema(bool overwrite=false) const;
};



class DBTable {
public:
    DBTable(const char* tablename, const char* primary, vector<string> fields, vector<string> types={});
    DBTable(DBConnection* database, const char* tablename, const char* primary, vector<string> fields, vector<string> types={});
    DBTable(DBConnection* database, const DBSchema* schema, const char* tablename, const char* primary, vector<string> fields, vector<string> types={});
    DBTable(DBConnection* database, const char* schemaname, const char* tablename, const char* primary, vector<string> fields, vector<string> types={});
    ~DBTable()=default;

protected:
    string schema, table, primary;
    vector<string> fields;
    vector<string> fieldtypes;      // stores type and constraint e.g. TEXT UNIQUE
    DBConnection* database=nullptr;

public:
    // getters
    string Name() const;
    string SchemaName() { return schema; }
    string TableName() { return table; }
    string Key() const { return primary; }
    string KeyType() const { return "SERIAL PRIMARY KEY"; }
    string Field(unsigned int i) const { return fields[i]; }        // throws error if i too big
    int Field(string& field) const;
    string Fieldtype(unsigned int i) const { return fieldtypes[i]; }
    string Fieldtype(string& field) const { return fieldtypes[Field(field)]; }
    vector<string> Fields() const { return fields; }
    vector<string> Fieldtypes() const { return fieldtypes; }
    pqxx::connection* Conn() { return database->Conn(); }
    
    // setters
    void Fieldtype(string& field, string& type);
    void connect(DBConnection* db) { database=db; }

    // member functions
    string query_create_table() const;
    string query_add_field(const string& fieldname, const string& type="") const;
    virtual string query_insert_into_table() const;     // not implemented in base class

    map<string,string> get_fields();        // get fields and fieldtypes from database
    void get_fieldtypes();    // get fieldtypes from database

    void add_field(const string& fieldname, const string& type="");     // does not update database!
    void create_schema(DBConnection* database, const char* schemaname);
    void create_schema(const char* schemaname);

};

// special case: a table with a primary key, a name field, and a field for a binary gpkg file
class GPKGTable : public DBTable {
public:
    GPKGTable(const char* tablename, const char* primary, const char* filename, const char* data,
              bool unique=false, bool with_timestamp=true, const char* timename="time");
    GPKGTable(DBConnection* database, const char* tablename, const char* primary, const char* filename, const char* data,
              bool unique=false, bool with_timestamp=true, const char* timename="time");
    GPKGTable(DBConnection* database, const DBSchema* schema, const char* tablename, const char* primary, const char* filename, const char* data,
              bool unique=false, bool with_timestamp=true, const char* timename="time");
    GPKGTable(DBConnection* database, const char* schemaname, const char* tablename, const char* primary, const char* filename, const char* data,
              bool unique=false, bool with_timestamp=true, const char* timename="time");
    ~GPKGTable()=default;

    // remove inherited member functions
    string Field(unsigned int)=delete;
    string Fieldtype(unsigned int)=delete;
    string Fieldtype(string)=delete;
    void Fields()=delete;
    void Fieldtypes()=delete;

public:     // new getters
    string NameField() const { return fields[0]; };
    string DataField() const { return fields[1]; };
    string TimeField() const {
        if (with_timestamp)     return fields[2];
        else { cerr << "No Timestamp field" << endl; return ""; }
    }
    string NameType() const { return fieldtypes[0]; };
    string DataType() const { return fieldtypes[1]; };

    string query_insert_into_table() const override;
    string query_get_from_table(const char* gpkg_name) const;

    bool is_unique();    // find out from database if name is unique
    void set_unique(bool value, bool modify_db=false);      // set GPKGTable.unique; if modify_db=true adjust table in database

    void prepare_db_table();
    bool upload_gpkg(const char* inputfile, const char* name_in_db, bool delete_local=false);       // returns false if failed
    bool download_gpkg(const char* name_in_db, const char* outputfile);

private:
    bool unique=false;
    bool with_timestamp=true;

};


#endif //FLUIDS_DBCONNECTION_H
