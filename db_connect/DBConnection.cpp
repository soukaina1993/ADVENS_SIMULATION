//
// Created by cornelia.blanke on 10.06.2025.
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include "DBConnection.h"


DBConnection::DBConnection(const char* dbname, const char* user, const char* password,
                           const char* host, int port) : dbname(dbname), user(user),password(password), host(host), port(port) {
    try {
        conn = new pqxx::connection(query_DB());
        if (!conn->is_open()) {
            cerr << "Cannot connect to database." << endl;
            exit(1);
        }
    }
    catch (const std::exception &e) {
        cerr << "Error : " << e.what() << endl;
        cerr << "Cannot connect to database." << endl;
        exit(1);
    }
}

// DBConnection member functions
string DBConnection::query_DB() const {
    stringstream query;
    query << "dbname=" << dbname << " user=" << user << " password=" << password << " host=" << host << " port=" << port;
    return query.str();
}



// Constructor
DBSchema::DBSchema(const char* schemaname) : name(schemaname) {}

// DBSchema member functions
string DBSchema::query_create_schema(const bool overwrite) const {
    stringstream query;
    if (overwrite)
        query << "DROP SCHEMA IF EXISTS " << SchemaName() << " CASCADE ; ";
    query << "CREATE SCHEMA IF NOT EXISTS " << SchemaName();
    return query.str();
}




// Constructors
DBTable::DBTable(const char* tablename, const char* primary, vector<string> fields, vector<string> types) : table(tablename),
        primary(primary), fields(std::move(fields)), fieldtypes(std::move(types)) {}

//DBTable::DBTable(DBConnection* database, const char* tablename, const char* primary, vector<string> fields, vector<string> types) : database(database),
//        name(tablename), primary(primary), fields(std::move(fields)), fieldtypes(std::move(types)) {}
DBTable::DBTable(DBConnection* database, const char* tablename, const char* primary,
                 vector<string> fields, vector<string> types) : DBTable(tablename, primary, std::move(fields), std::move(types)) {
    this->database=database;
}

DBTable::DBTable(DBConnection* database, const DBSchema* schema, const char* tablename, const char* primary,
                 vector<string> fields, vector<string> types) : DBTable(database, tablename,
                                                                        primary, std::move(fields), std::move(types)) {
    this->schema = schema->SchemaName();
}

DBTable::DBTable(DBConnection* database, const char* schemaname, const char* tablename, const char* primary,
                 vector<string> fields, vector<string> types) : DBTable(database, tablename,
                                                                        primary, std::move(fields), std::move(types)) {
    create_schema(schemaname);
    schema = schemaname;
}





GPKGTable::GPKGTable(const char* tablename, const char* primary, const char* filename, const char* data,
                     bool unique, bool with_timestamp, const char* timename) :
                         DBTable(tablename, primary, {filename, data}, {"TEXT UNIQUE", "BYTEA"}),
                         unique(unique),
                         with_timestamp(with_timestamp) {
    if (!this->unique)
        fieldtypes[0] = "TEXT";
    if (this->with_timestamp)
        add_field(timename, "TIMESTAMP");
}

GPKGTable::GPKGTable(DBConnection* database, const char* tablename, const char* primary, const char* filename, const char* data,
                     bool unique, bool with_timestamp, const char* timename) :
                        DBTable(database, tablename, primary, {filename, data}, {"TEXT UNIQUE", "BYTEA"}),
                        unique(unique),
                        with_timestamp(with_timestamp) {
    if (!this->unique)
        fieldtypes[0] = "TEXT";
    if (this->with_timestamp)
        add_field(timename, "TIMESTAMP");
}

GPKGTable::GPKGTable(DBConnection* database, const DBSchema* schema, const char* tablename, const char* primary, const char* filename, const char* data,
                    bool unique, bool with_timestamp, const char* timename) :
                        DBTable(database, schema, tablename, primary, {filename, data}, {"TEXT UNIQUE", "BYTEA"}),
                        unique(unique),
                        with_timestamp(with_timestamp) {
    if (!this->unique)
        fieldtypes[0] = "TEXT";
    if (this->with_timestamp)
        add_field(timename, "TIMESTAMP");
}

GPKGTable::GPKGTable(DBConnection* database, const char* schemaname, const char* tablename, const char* primary, const char* filename, const char* data,
                    bool unique, bool with_timestamp, const char* timename) :
                        DBTable(database, schemaname, tablename, primary, {filename, data}, {"TEXT UNIQUE", "BYTEA"}),
                        unique(unique),
                        with_timestamp(with_timestamp) {
    if (!this->unique)
        fieldtypes[0] = "TEXT";
    if (this->with_timestamp)
        add_field(timename, "TIMESTAMP");
}



// DBTable member functions
string DBTable::Name() const {
    if (schema.empty())
        return table;
    else
        return schema + "." + table;
}

int DBTable::Field(string& field) const {
    auto it = find(fields.begin(), fields.end(), field);
    if (it != fields.end())
        return (int)(it - fields.begin());
    else {
        cerr << "Field " << field << " was not found" << endl;
        exit(1);
    }
}

void DBTable::Fieldtype(string& field, string& type) {
    if (fieldtypes.size() != fields.size())
        fieldtypes.resize(fields.size());
    fieldtypes[Field(field)] = type;
}

string DBTable::query_insert_into_table() const {
    cerr << "Not implemented" << endl;
    exit(1);
}

string DBTable::query_create_table() const {
    if (fields.size() != fieldtypes.size()){
         cerr << "Cannot create table: Specify correct number of fields with their field types" << endl;
         exit(1);
    }
    stringstream query;
    query << "CREATE TABLE IF NOT EXISTS " << Name() << " (";
    query << Key()  << " " << KeyType();
    for (int i=0; i<fields.size(); i++)
        query << ", " << fields[i] << " " << fieldtypes[i];
    query << ")";
    return query.str();
}

string DBTable::query_add_field(const string& fieldname, const string& type) const {
    stringstream query;
    query << "ALTER TABLE " << Name() << " ADD " << fieldname << " " << type;
    return query.str();
}

map<string, string> DBTable::get_fields() {
    if (database == nullptr){
        cerr << "No database specified" << endl;
        exit(1);
    }
    map<string, string> dict;
    try {
        stringstream query;
        query << "SELECT column_name, data_type FROM information_schema.columns WHERE table_name = '";
        query << Name() << "'";
        //query << " AND column_name = '" << field << "'";
        pqxx::work txn(*Conn());
        pqxx::result r = txn.exec(query.str());

        for (auto && row : r) {
            if (strcasecmp(primary.c_str(), row[0].as<std::string>().c_str()) == 0)
                continue;
            auto field = row[0].as<std::string>();
            auto type = row[1].as<std::string>();
            dict[field] = type;     // add entry to map
        }
    }
    catch (const std::exception &e) {
        cerr << "Error : " << e.what() << endl;
        cerr << "Fields cannot be fetched" << endl;
    }
    return dict;
}

void DBTable::get_fieldtypes() {
    map<string, string> dict = get_fields();
    for (auto& field : fields)
        Fieldtype(field, dict.at(field));      // set fieldtype (error if field is not found)
}

void DBTable::add_field(const string& fieldname, const string& type){
    auto it=find(fields.begin(), fields.end(), fieldname);
    if (it == fields.end()) {      // not yet there
        fields.push_back(fieldname);
        if (!type.empty()) {
            fieldtypes.resize(fields.size());   // works in any case
            fieldtypes.back() = type;
        }
    }
    else {
        ptrdiff_t pos = it - fields.begin();
        if (fieldtypes.size() > pos && fieldtypes[pos] != type)
            cerr << "Field " << fieldname << "is already defined as " << fieldtypes[pos] << endl;
    }
}

void DBTable::create_schema(DBConnection* db, const char* schemaname){
    if (db != nullptr) {
        DBSchema schema(schemaname);
        pqxx::work txn(*db->Conn());
        txn.exec(schema.query_create_schema());
        txn.commit();
    }
    else {
        cerr << "No database connected. Cannot create schema " << schemaname << endl;
        exit(1);
    }
}

void DBTable::create_schema(const char* schemaname){
    create_schema(database, schemaname);
}

// GPKGTable member functions
string GPKGTable::query_insert_into_table() const {
    stringstream query;
    query << "INSERT INTO " << Name();
    if (!with_timestamp)
        query << " (" << NameField() <<", " << DataField() << ") VALUES ($1, $2)";
    else
        query << " (" << NameField() <<", " << DataField() << ", " << TimeField() << ") VALUES ($1, $2, NOW())";
    if (unique)     // add ON CONFLICT
        query << " ON CONFLICT (" << NameField() << ") DO UPDATE SET " << DataField() << " = EXCLUDED." << DataField();      // update if unique = true
    return query.str();
}

string GPKGTable::query_get_from_table(const char* gpkg_name) const {
    stringstream query;
    query << "SELECT " << NameField() << ", " << DataField();
    query << " FROM " << Name();
    query << " WHERE " << NameField() << " = '" << gpkg_name << "'";
    // get latest version by timestamp or key
    query << " ORDER BY " << (with_timestamp ? TimeField() : Key()) << " DESC LIMIT 1";
    return query.str();
}

bool GPKGTable::is_unique() {
    try {
        stringstream query;
        query << "SELECT COUNT(*) FROM INFORMATION_SCHEMA.TABLE_CONSTRAINTS WHERE constraint_name = '";
        query << Name() << "_" << NameField() << "_key' AND constraint_type = 'UNIQUE'";

        pqxx::work txn(*Conn());
        pqxx::result r = txn.exec(query.str());
        auto value = r[0][0].as<int>();
        if (value > 0)
            return true;
    }
    catch (const std::exception &e) {
        cerr << "Error : " << e.what() << endl;
        cerr << "Uniqueness cannot be fetched" << endl;
        return false;
    }
    return false;
}

void GPKGTable::set_unique(const bool value, const bool modify_db){
    unique = value;
    unique ? NameType()="TEXT UNIQUE" : NameType()="TEXT";

    if (modify_db){
        if (Conn() == nullptr)
            cerr << "No connected database" << endl;
        else {
            stringstream query;
            if (value)
                query << "ALTER TABLE " << Name() << " ADD CONSTRAINT " << Name() << "_" << NameField() << "_key UNIQUE (" << NameField() << ")";
            else if (is_unique())   // is unique on database but should not be
                query << "ALTER TABLE " << Name() << " DROP CONSTRAINT " << Name() << "_" << NameField() << "_key";
            pqxx::nontransaction txn(*Conn());
            txn.exec(query.str());
        }
    }
}

void GPKGTable::prepare_db_table() {
    pqxx::nontransaction txn(*Conn());
    stringstream query;
    query << "SELECT COUNT(*) FROM pg_tables WHERE  tablename = LOWER('" << TableName() << "')";
    if (!schema.empty())
        query << " AND schemaname  = LOWER('" << SchemaName() << "')";
    pqxx::result r = txn.exec(query.str());
    bool table_exists = r[0][0].as<bool>();
    if (!table_exists)
        txn.exec(query_create_table());
    else {    // check for missing fields
        stringstream().swap(query);     // reset and clear query
        query << "SELECT column_name, data_type FROM information_schema.columns WHERE table_name = LOWER('";
        query << TableName() << "')";
        if (!schema.empty())
            query << "AND table_schema = LOWER('" << SchemaName() << "')";
        r = txn.exec(query.str());
        for (int i=0; i<fields.size(); i++){    // loop over name of fields
            bool found = false;
            for (auto && row : r) {
                if (strcasecmp(fields[i].c_str(), row[0].as<std::string>().c_str()) == 0) {
                    found = true;
                    // check if fieldtype is correct
                    auto type_in_db = row[1].as<std::string>();
                    int minlen = min((int) fieldtypes[i].size(), (int) type_in_db.size());
                    string fieldtype = fieldtypes[i].substr(0, minlen);
                    type_in_db = type_in_db.substr(0, minlen);

                    if(strcasecmp(fieldtype.c_str(), type_in_db.c_str()) != 0) {
                        cerr << "Field " << fields[i] << " has wrong type " << fieldtype << " " << type_in_db << endl;
                        exit(1);
                    }
                    else
                        break;  // check next field
                }
            }
            if (!found)       // create new field on table
                txn.exec(query_add_field(fields[i], fieldtypes[i]));
        }
    }
    txn.commit();
    set_unique(unique, true);       // add or remove UNIQUE constraint
}

bool GPKGTable::upload_gpkg(const char* inputfile, const char* name_in_db, const bool delete_local) {
    try {
        string fileData;
        {
            ifstream file(inputfile, std::ios::binary);
            ostringstream oss;
            oss << file.rdbuf();
            fileData = oss.str();
        }   // with {..} the local variables will be deleted and the inputfile will be closed
        prepare_db_table();

        pqxx::work txn(*Conn());
        txn.exec(query_insert_into_table(),
                 pqxx::params(name_in_db, pqxx::binary_cast(fileData))
        );
        txn.commit();
    }
    catch (const std::exception &e) {
        cerr << "Error : " << e.what() << endl;
        cerr << inputfile << " was not uploaded" << endl;
        return false;
    }

    if (delete_local)
        std::remove(inputfile);

    return true;
}

bool GPKGTable::download_gpkg(const char* name_in_db, const char* outputfile) {
    try {
        // start transaction
        pqxx::work txn(*Conn());
        pqxx::result r = txn.exec(query_get_from_table(name_in_db));

        if (r.empty()) {
            cerr << name_in_db << " was not found in database." << endl;
            return false;
        }

        auto binary_data = r[0]["data"].as<std::basic_string<std::byte>>();

        std::ofstream file(outputfile, std::ios::binary);
        for(auto s : binary_data)
            file << (char)s;
        file.close();
    }
    catch (const std::exception &e) {
        cerr << "Error : " << e.what() << endl;
        cerr << name_in_db << " was not downloaded" << endl;
        return false;
    }
    return true;
}
